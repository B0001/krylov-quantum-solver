"""Gate for provenance-based dataset compilation (bead chem-fdq).

compile_from_gcs.py used to recover each result's geometry by computing local
CASCI energies for every candidate coordinate and taking the argmin. That made
the coordinate label a function of the reference, so the compiled dataset
agreed with CASCI to <1e-4 Ha by construction. These tests pin the replacement:
the producer states what it ran, the consumer reads it, and the CASCI
comparison is then free to disagree.
"""
import importlib
import json
import sys

import pytest

fastapi_testclient = pytest.importorskip("fastapi.testclient")
pytest.importorskip("google.cloud.storage")
TestClient = fastapi_testclient.TestClient

import compile_from_gcs as cfg  # noqa: E402 -- must follow the importorskip guards

KEY = "provenance-test-key"
BODY = {
    "molecule": "H 0 0 0; H 0 0 0.74",
    "basis": "sto-3g",
    "active_space": [2, 2],
    "mode": "certified",
    "krylov_dim": 10,
}


class RecordingRedis:
    """Captures what the API writes, so the job record can be inspected."""

    def __init__(self):
        self.store = {}
        self.queue = []
        self.counts = {}

    def setex(self, key, ttl, val):
        self.store[key] = val

    def set(self, key, val):
        self.store[key] = val

    def get(self, key):
        return self.store.get(key)

    def lpush(self, key, val):
        self.queue.append(val)

    def ping(self):
        return True

    def pipeline(self):
        return _Pipe(self)


class _Pipe:
    def __init__(self, parent):
        self.parent, self.keys = parent, []

    def incr(self, k):
        self.keys.append(k)

    def expire(self, k, ttl):
        pass

    def execute(self):
        out = []
        for k in self.keys:
            self.parent.counts[k] = self.parent.counts.get(k, 0) + 1
            out.append(self.parent.counts[k])
        return out


@pytest.fixture
def api(monkeypatch):
    monkeypatch.setenv("ALLOWED_API_KEYS", KEY)
    monkeypatch.delenv("RESULTS_BUCKET", raising=False)
    sys.modules.pop("src.certchem.api", None)
    mod = importlib.import_module("src.certchem.api")
    mod.redis_client = RecordingRedis()
    return mod


# --- the API must carry caller labels through to the job ---------------------

def test_metadata_is_persisted_on_the_job(api):
    c = TestClient(api.app)
    meta = {"family": "n2_curve", "coordinate": 1.4}
    r = c.post("/v1/energy", json={**BODY, "metadata": meta}, headers={"X-API-Key": KEY})
    assert r.status_code == 202, r.text

    job = json.loads(api.redis_client.store[f"job:{r.json()['job_id']}"])
    assert job["request"]["metadata"] == meta
    assert job["request"]["molecule"] == BODY["molecule"]


def test_metadata_is_optional(api):
    c = TestClient(api.app)
    r = c.post("/v1/energy", json=BODY, headers={"X-API-Key": KEY})
    assert r.status_code == 202
    job = json.loads(api.redis_client.store[f"job:{r.json()['job_id']}"])
    assert job["request"]["metadata"] is None


def test_oversized_metadata_is_rejected_upfront(api):
    c = TestClient(api.app)
    fat = {"junk": "x" * (api.METADATA_MAX_BYTES + 100)}
    r = c.post("/v1/energy", json={**BODY, "metadata": fat}, headers={"X-API-Key": KEY})
    assert r.status_code == 422
    assert api.redis_client.queue == [], "a rejected job must not reach the queue"


# --- the compile path must read provenance, never infer it -------------------

def blob(**over):
    data = {
        "provenance": {
            "job_id": "abc123",
            "molecule": "N 0 0 0; N 0 0 1.4",
            "basis": "6-31g",
            "active_space": [6, 6],
            "mode": "certified",
            "krylov_dim": 10,
            "metadata": {"family": "n2_curve", "coordinate": 1.4},
        },
        "best_estimate_hartree": -109.1,
        "lower_bound_hartree": -109.2,
        "upper_bound_hartree": -109.0,
        "bracket_width_hartree": 0.2,
        "wall_time_s": 3.5,
    }
    data.update(over)
    return data


def test_row_is_built_from_recorded_provenance():
    row = cfg.row_from_blob(blob())
    assert row["family"] == "n2_curve"
    assert row["coordinate_var"] == 1.4
    assert row["molecule"] == "N 0 0 0; N 0 0 1.4"
    assert row["basis"] == "6-31g"
    assert row["active_space"] == "6,6"
    assert row["job_id"] == "abc123"


def test_blob_without_provenance_is_refused_not_guessed():
    """The old code would have reverse-engineered a label from the energy."""
    legacy = {"best_estimate_hartree": -109.1}
    assert cfg.row_from_blob(legacy) is None
    assert cfg.row_from_blob(blob(provenance={})) is None


def test_missing_metadata_still_keeps_the_geometry():
    """No caller label is not the same as no provenance -- keep what we know."""
    b = blob()
    b["provenance"]["metadata"] = None
    row = cfg.row_from_blob(b)
    assert row["family"] == "unlabelled"
    assert row["coordinate_var"] is None
    assert row["molecule"] == "N 0 0 0; N 0 0 1.4"


@pytest.mark.parametrize(
    "gone",
    ["identify_family_by_energy", "get_standard_geometries", "calculate_local_reference_energy"],
)
def test_energy_based_identification_is_gone(gone):
    """The circular helpers must not come back.

    Asserted on the module's attributes rather than by grepping its source:
    the docstring legitimately names what was removed and why, and a grep
    cannot tell an explanation from a reintroduction.
    """
    assert not hasattr(cfg, gone)


# --- the CASCI comparison must now be able to disagree -----------------------

def test_verify_computes_a_real_reference_and_accepts_a_true_bracket():
    pytest.importorskip("pyscf")
    row = cfg.row_from_blob(blob(
        provenance={
            "job_id": "h2", "molecule": "H 0 0 0; H 0 0 0.74", "basis": "sto-3g",
            "active_space": [2, 2], "mode": "certified", "krylov_dim": 10,
            "metadata": {"family": "h2", "coordinate": 0.74},
        },
        best_estimate_hartree=-1.1373, lower_bound_hartree=-1.2,
        upper_bound_hartree=-1.0, bracket_width_hartree=0.2,
    ))
    out = cfg.verify(row)
    assert out["casci_reference_hartree"] == pytest.approx(-1.1373, abs=1e-3)
    assert out["bracket_contains_reference"] is True


def test_verify_can_fail_a_bracket_that_excludes_the_reference():
    """The check earns its name only if it can come back False."""
    pytest.importorskip("pyscf")
    row = cfg.row_from_blob(blob(
        provenance={
            "job_id": "h2", "molecule": "H 0 0 0; H 0 0 0.74", "basis": "sto-3g",
            "active_space": [2, 2], "mode": "certified", "krylov_dim": 10,
            "metadata": {"family": "h2", "coordinate": 0.74},
        },
        best_estimate_hartree=-0.5, lower_bound_hartree=-0.6,
        upper_bound_hartree=-0.4, bracket_width_hartree=0.2,
    ))
    out = cfg.verify(row)
    assert out["bracket_contains_reference"] is False
    assert out["error_vs_casci_hartree"] > 0.5


# --- the worker is what actually stamps the blob -----------------------------

def load_worker(monkeypatch):
    monkeypatch.delenv("RESULTS_BUCKET", raising=False)
    sys.modules.pop("src.certchem.worker", None)
    return importlib.import_module("src.certchem.worker")


def run_one_job(monkeypatch, request, solver_result):
    """Drive process_job with a stubbed solver; return the uploaded blob."""
    w = load_worker(monkeypatch)
    w.redis_client = RecordingRedis()
    w.redis_client.store["job:J1"] = json.dumps({"job_id": "J1", "status": "queued",
                                                 "request": request})
    uploaded = {}
    monkeypatch.setattr(w, "certified_energy", lambda **kw: solver_result)
    monkeypatch.setattr(w, "upload_result_to_gcs",
                        lambda job_id, data: uploaded.update(data))
    w.process_job("J1")
    assert json.loads(w.redis_client.store["job:J1"])["status"] == "completed"
    return uploaded


REQUEST = {
    "molecule": "N 0 0 0; N 0 0 1.4",
    "basis": "6-31g",
    "active_space": [6, 6],
    "mode": "fast",
    "krylov_dim": 10,
    "metadata": {"family": "n2_curve", "coordinate": 1.4},
}


def test_worker_stamps_provenance_onto_the_result_blob(monkeypatch):
    uploaded = run_one_job(monkeypatch, REQUEST, -109.1)
    prov = uploaded["provenance"]
    assert prov["molecule"] == REQUEST["molecule"]
    assert prov["basis"] == "6-31g"
    assert prov["active_space"] == [6, 6]
    assert prov["metadata"] == {"family": "n2_curve", "coordinate": 1.4}
    assert prov["job_id"] == "J1"


def test_stamped_blob_round_trips_through_the_compiler(monkeypatch):
    """End to end: what the worker writes is what the compiler can label."""
    uploaded = run_one_job(monkeypatch, REQUEST, -109.1)
    row = cfg.row_from_blob(json.loads(json.dumps(uploaded)))
    assert row is not None
    assert row["family"] == "n2_curve"
    assert row["coordinate_var"] == 1.4
    assert row["molecule"] == REQUEST["molecule"]
    assert row["wall_time_s"] is not None


def test_worker_handles_a_job_submitted_without_metadata(monkeypatch):
    uploaded = run_one_job(monkeypatch, {**REQUEST, "metadata": None}, -109.1)
    assert uploaded["provenance"]["metadata"] is None
    assert uploaded["provenance"]["molecule"] == REQUEST["molecule"]
