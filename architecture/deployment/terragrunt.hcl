# Terragrunt configuration for CertChem infrastructure deployment
#
# Enables dry backend management, dynamic Google provider configuration,
# and keeps deployment environment inputs separated.
#
# Usage:
#   export GCP_PROJECT_ID="your-gcp-project-id"
#   export GCP_IMAGE_URL="us-docker.pkg.dev/your-gcp-project-id/certchem/app:latest"
#   terragrunt apply

# 1. Source local Terraform module (current directory)
terraform {
  source = "."
}

# 2. Automatically generate the backend.tf file for GCS remote state
remote_state {
  backend = "gcs"
  generate = {
    path      = "backend.tf"
    if_exists = "overwrite_terragrunt"
  }
  config = {
    project  = get_env("GCP_PROJECT_ID", "unconfigured-project-id")
    bucket   = "${get_env("GCP_PROJECT_ID", "unconfigured-project-id")}-certchem-tf-state"
    prefix   = "state"
    location = get_env("GCP_REGION", "us-central1")
  }
}

# 3. Supply inputs to the underlying Terraform variables declared in main.tf
inputs = {
  project_id = get_env("GCP_PROJECT_ID", "unconfigured-project-id")
  region     = get_env("GCP_REGION", "us-central1")
  image      = get_env("GCP_IMAGE_URL", "us-docker.pkg.dev/unconfigured-project-id/certchem/app:latest")
}
