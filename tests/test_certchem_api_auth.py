"""Gate for the ADR-0009 app-layer auth: API keys + per-key rate limits.

Not an SDD spec gate (no scientific claim, no reference energy), so it is
deliberately NOT named test_*_spec.py -- that glob is the physics gates the
Makefile runs in isolated processes.

The API module reads ALLOWED_API_KEYS and builds its Redis client at import
time, so every test imports it fresh under a patched environment with a fake
Redis standing in for Memorystore.
"""
import importlib
import sys

import pytest

fastapi_testclient = pytest.importorskip("fastapi.testclient")
TestClient = fastapi_testclient.TestClient

GOOD_KEY = "test-key-alpha"


class FakeRedis:
    """Minimal INCR/EXPIRE/pipeline stand-in; counters only, no TTL clock."""

    def __init__(self):
        self.counts = {}

    def pipeline(self):
        return FakePipeline(self)

    # submit_energy_job touches these; they only need to not explode.
    def setex(self, *a, **k):
        return True

    def lpush(self, *a, **k):
        return 1

    def get(self, *a, **k):
        return None

    def ping(self):
        return True


class FakePipeline:
    def __init__(self, parent):
        self.parent = parent
        self.queued = []

    def incr(self, key):
        self.queued.append(key)

    def expire(self, key, ttl):
        pass

    def execute(self):
        out = []
        for key in self.queued:
            self.parent.counts[key] = self.parent.counts.get(key, 0) + 1
            out.append(self.parent.counts[key])
        return out


def load_api(monkeypatch, allowed=GOOD_KEY, limit="100"):
    monkeypatch.setenv("ALLOWED_API_KEYS", allowed)
    monkeypatch.setenv("RATE_LIMIT_PER_MIN", limit)
    monkeypatch.delenv("RESULTS_BUCKET", raising=False)
    sys.modules.pop("src.certchem.api", None)
    api = importlib.import_module("src.certchem.api")
    api.redis_client = FakeRedis()
    return api


@pytest.fixture
def client(monkeypatch):
    return TestClient(load_api(monkeypatch).app)


ENERGY_BODY = {
    "molecule": "H 0 0 0; H 0 0 0.74",
    "basis": "sto-3g",
    "active_space": [2, 2],
    "mode": "certified",
    "krylov_dim": 10,
}


def test_public_endpoints_need_no_key(client):
    """GCP health checks and the limits doc must work without credentials."""
    assert client.get("/health").status_code == 200
    assert client.get("/v1/limits").status_code == 200


@pytest.mark.parametrize("headers", [{}, {"X-API-Key": "wrong-key"}])
def test_job_endpoints_reject_missing_or_bad_key(client, headers):
    assert client.post("/v1/energy", json=ENERGY_BODY, headers=headers).status_code == 401
    assert client.get("/v1/jobs/deadbeef", headers=headers).status_code == 401


def test_good_key_is_admitted(client):
    """A valid key gets past auth: 202 queued, not 401."""
    r = client.post("/v1/energy", json=ENERGY_BODY, headers={"X-API-Key": GOOD_KEY})
    assert r.status_code == 202, r.text
    assert r.json()["status"] == "queued"


def test_rate_limit_trips_after_quota(monkeypatch):
    """The Nth+1 request in a window is refused with 429 and a Retry-After."""
    api = load_api(monkeypatch, limit="3")
    c = TestClient(api.app)
    h = {"X-API-Key": GOOD_KEY}
    for i in range(3):
        assert c.get(f"/v1/jobs/job{i}", headers=h).status_code == 404  # past auth
    r = c.get("/v1/jobs/job4", headers=h)
    assert r.status_code == 429
    assert "Retry-After" in r.headers


def test_rate_limit_is_per_key(monkeypatch):
    """One key burning its quota must not throttle a different key."""
    api = load_api(monkeypatch, allowed=f"{GOOD_KEY},second-key", limit="2")
    c = TestClient(api.app)
    for i in range(2):
        c.get(f"/v1/jobs/j{i}", headers={"X-API-Key": GOOD_KEY})
    assert c.get("/v1/jobs/x", headers={"X-API-Key": GOOD_KEY}).status_code == 429
    assert c.get("/v1/jobs/x", headers={"X-API-Key": "second-key"}).status_code == 404


def test_unconfigured_deployment_fails_closed(monkeypatch):
    """No ALLOWED_API_KEYS must serve 503, never an open unauthenticated API."""
    api = load_api(monkeypatch, allowed="")
    c = TestClient(api.app)
    r = c.post("/v1/energy", json=ENERGY_BODY, headers={"X-API-Key": GOOD_KEY})
    assert r.status_code == 503
    assert c.post("/v1/energy", json=ENERGY_BODY).status_code == 503


def test_raw_keys_are_not_stored_in_redis(monkeypatch):
    """Rate-limit buckets are keyed by digest, so Redis never holds a raw key."""
    api = load_api(monkeypatch)
    c = TestClient(api.app)
    c.get("/v1/jobs/j", headers={"X-API-Key": GOOD_KEY})
    buckets = list(api.redis_client.counts)
    assert buckets, "expected a rate-limit bucket to be created"
    assert all(GOOD_KEY not in b for b in buckets), buckets
