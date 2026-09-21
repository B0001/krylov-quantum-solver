# CertChem infrastructure — Terraform for v1 (ADR-0009)
terraform {
  required_providers {
    google = {
      source  = "hashicorp/google"
      version = "~> 6.0"
    }
  }
}

variable "project_id" {
  type = string
}

variable "region" {
  type    = string
  default = "us-central1"
}

variable "image" {
  type = string # e.g. us-docker.pkg.dev/PROJECT/certchem/app:GIT_SHA
}

provider "google" {
  project = var.project_id
  region  = var.region
}

# --- State stores -----------------------------------------------------------
resource "google_storage_bucket" "results" {
  name                        = "${var.project_id}-certchem-results"
  location                    = var.region
  uniform_bucket_level_access = true
  versioning {
    enabled = true
  }
}

resource "google_redis_instance" "queue" {
  name           = "certchem-redis"
  tier           = "BASIC" # queue+cache are reconstructible; no HA needed for v1
  memory_size_gb = 1
  region         = var.region
}

resource "google_vpc_access_connector" "conn" {
  name          = "certchem-vpc-conn"
  region        = var.region
  network       = "default"
  ip_cidr_range = "10.8.0.0/28"
  min_instances = 2
  max_instances = 3
}

# --- Identities (least privilege) --------------------------------------------
resource "google_service_account" "api" {
  account_id   = "certchem-api-sa"
  display_name = "CertChem API Service Account"
}

resource "google_service_account" "worker" {
  account_id   = "certchem-worker-sa"
  display_name = "CertChem Worker Service Account"
}

resource "google_storage_bucket_iam_member" "api_read" {
  bucket = google_storage_bucket.results.name
  role   = "roles/storage.objectViewer"
  member = "serviceAccount:${google_service_account.api.email}"
}

resource "google_storage_bucket_iam_member" "worker_write" {
  bucket = google_storage_bucket.results.name
  role   = "roles/storage.objectAdmin"
  member = "serviceAccount:${google_service_account.worker.email}"
}

# --- API-key material (ADR-0009: app-layer auth) -------------------------------
# The secret RESOURCE is declared here; its VALUE is never in Terraform, so the
# keys stay out of tfstate. Populate out of band:
#   printf 'key1,key2' | gcloud secrets versions add certchem-api-keys --data-file=-
resource "google_secret_manager_secret" "api_keys" {
  secret_id = "certchem-api-keys"
  replication {
    auto {}
  }
}

resource "google_secret_manager_secret_iam_member" "api_keys_access" {
  secret_id = google_secret_manager_secret.api_keys.id
  role      = "roles/secretmanager.secretAccessor"
  member    = "serviceAccount:${google_service_account.api.email}"
}

# --- API service --------------------------------------------------------------
resource "google_cloud_run_v2_service" "api" {
  name                = "certchem-api"
  location            = var.region
  deletion_protection = false

  template {
    service_account = google_service_account.api.email

    scaling {
      min_instance_count = 0
      max_instance_count = 10
    }

    max_instance_request_concurrency = 20

    vpc_access {
      connector = google_vpc_access_connector.conn.id
      egress    = "PRIVATE_RANGES_ONLY" # Redis only; API has no other egress needs
    }

    containers {
      image = var.image

      # $PORT is injected by Cloud Run; app binds 0.0.0.0:$PORT (ADR-0009)
      env {
        name  = "REDIS_HOST"
        value = google_redis_instance.queue.host
      }
      env {
        name  = "RESULTS_BUCKET"
        value = google_storage_bucket.results.name
      }
      env {
        name  = "SYNC_TIMEOUT_S"
        value = "600"
      }
      # Comma-separated allowed keys. Absent => the API fails closed with 503
      # rather than serving unauthenticated; see require_api_key in api.py.
      env {
        name = "ALLOWED_API_KEYS"
        value_source {
          secret_key_ref {
            secret  = google_secret_manager_secret.api_keys.secret_id
            version = "latest"
          }
        }
      }
      env {
        name  = "RATE_LIMIT_PER_MIN"
        value = "100"
      }

      resources {
        limits = {
          cpu    = "1"
          memory = "1Gi"
        }
      }
    }

    timeout = "600s"
  }
}

# Public HTTPS ingress; auth is application-layer API keys + rate limits
resource "google_cloud_run_v2_service_iam_member" "public" {
  name     = google_cloud_run_v2_service.api.name
  location = var.region
  role     = "roles/run.invoker"
  member   = "allUsers"
}

# --- Worker (queue consumer, no ingress) ---------------------------------------
resource "google_cloud_run_v2_job" "worker" {
  name                = "certchem-worker"
  location            = var.region
  deletion_protection = false

  template {
    template {
      service_account = google_service_account.worker.email

      vpc_access {
        connector = google_vpc_access_connector.conn.id
        egress    = "PRIVATE_RANGES_ONLY" # chemistry needs no internet
      }

      containers {
        image   = var.image
        command = ["python", "-m", "src.certchem.worker"]

        env {
          name  = "REDIS_HOST"
          value = google_redis_instance.queue.host
        }
        env {
          name  = "RESULTS_BUCKET"
          value = google_storage_bucket.results.name
        }
        env {
          name  = "QUEUED_TIMEOUT_S"
          value = "7200"
        }

        resources {
          limits = {
            cpu    = "4"
            memory = "16Gi"
          }
        }
      }

      max_retries = 1
      timeout     = "7200s"
    }

    parallelism = 5
  }
}

output "api_url" {
  value = google_cloud_run_v2_service.api.uri
}

output "redis_host" {
  value = google_redis_instance.queue.host
}
