"""
scripts/run_gates.sh must distinguish a KILLED gate from a FAILED one.

Before issue #36 every non-zero exit printed the same bare "FAIL", so a block2 SIGSEGV -- the
exact hazard the per-process gate isolation exists to prevent (CLAUDE.md) -- was indistinguishable
from an ordinary assertion failure, and a process killed before writing left a 0-byte log the FAIL
line pointed at.

Not named *_spec.py on purpose: it is infrastructure, not a physics claim, and the `test_*_spec.py`
glob is what run_gates.sh itself iterates -- a gate that runs the gate runner would recurse.
Deliberately signal-based (os.kill) rather than a real crash: the point is the classification, and
a genuine null deref is not portable enough to gate on.
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest

RUNNER = Path(__file__).resolve().parent.parent / "scripts" / "run_gates.sh"

GATES = {
    "test_green.py": "def test_ok():\n    assert True\n",
    "test_assertfail.py": "def test_bad():\n    assert False\n",
    "test_segv.py": "import os, signal\n\nos.kill(os.getpid(), signal.SIGSEGV)\n",
    "test_oom.py": "import os, signal\n\nos.kill(os.getpid(), signal.SIGKILL)\n",
}


@pytest.fixture(scope="module")
def report(tmp_path_factory):
    """Run the real runner over synthetic gates in a scratch cwd, return its stdout."""
    work = tmp_path_factory.mktemp("gaterunner")
    for name, src in GATES.items():
        (work / name).write_text(src)
    env = {**os.environ, "GATE_RUN": "", "GATE_JOBS": "1",
           "GATE_GLOB": f"{work}/test_*.py", "GATE_NO_CACHE": "1"}
    # `$RUN python -m pytest` with RUN="" needs `python` to be this interpreter.
    env["PATH"] = f"{Path(sys.executable).parent}{os.pathsep}{env['PATH']}"
    proc = subprocess.run(["bash", str(RUNNER)], cwd=work, env=env,
                          capture_output=True, text=True, timeout=300)
    return {ln.split()[1].rsplit("/", 1)[-1]: ln
            for ln in proc.stdout.splitlines() if ln.startswith(("PASS", "FAIL"))}


def test_passing_gate_still_passes(report):
    assert report["test_green.py"].startswith("PASS")


def test_assertion_failure_reports_pytest_exit_code(report):
    """A real test failure is exit 1 -- and must NOT be labelled as a signal."""
    line = report["test_assertfail.py"]
    assert line.startswith("FAIL") and "exit=1" in line, line
    assert "SIG" not in line, line


def test_segfault_is_named_not_silently_failed(report):
    """THE POINT: 139 must say SIGSEGV and point at the block2 isolation, not read as a failure."""
    line = report["test_segv.py"]
    assert "exit=139" in line and "SIGSEGV" in line, line
    assert "block2" in line, line


def test_oom_kill_is_named(report):
    line = report["test_oom.py"]
    assert "exit=137" in line and "SIGKILL" in line, line


def test_empty_log_is_called_out(report):
    """Killed before writing a byte -> say so, or FAIL points at an empty file."""
    assert "died before writing" in report["test_oom.py"], report["test_oom.py"]
    # an assertion failure DOES write output, so it must not claim an empty log
    assert "died before writing" not in report["test_assertfail.py"]
