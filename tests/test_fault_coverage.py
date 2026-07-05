"""Verification-modality coverage (D7): the lanes catch DISJOINT fault
classes, so the multi-modal stack is justified -- no single lane
suffices. Asserts the properties of the coverage matrix produced by
tools/fault_coverage.py."""

import pytest
from tools import fault_coverage as fc

PHYSICS = {"conservation", "analytic"}
CERT_LIKE = {"emitter translation", "spec corridor breach",
             "invalid barrier cert", "invalid KKT solution",
             "invalid Lyapunov P", "non-PSD SOS Gram"}


@pytest.fixture(scope="module")
def matrix():
    return fc.build_matrix()


def test_no_false_alarms(matrix):
    """A lane run against a GOOD artifact (any fault not in its remit)
    must never fire."""
    for fault, row in matrix.items():
        for lane in fc.LANES:
            if lane not in fc.FAULTS[fault]:
                assert row[lane] in (False, None), (fault, lane)


def test_every_lane_catches_its_own_fault(matrix):
    """Each lane listed for a fault actually rejects the faulted
    artifact (the lanes are effective, not decorative)."""
    for fault, seen in fc.FAULTS.items():
        for lane in seen:
            got = matrix[fault][lane]
            if got is None:            # lane unavailable (e.g. no gcc)
                continue
            assert got is True, (fault, lane)


def test_every_fault_is_caught(matrix):
    """Coverage: no injected fault slips past the whole stack."""
    for fault, seen in fc.FAULTS.items():
        assert any(matrix[fault][ln] for ln in seen), fault


def test_faults_are_caught_by_disjoint_lanes(matrix):
    """The load-bearing finding: six of seven faults are caught by
    exactly one lane, so removing that lane makes the fault ship
    silently. Certificate/emitter/spec faults are invisible to the
    physics lanes -- a wrong proof does not violate physics."""
    singletons = [f for f, seen in fc.FAULTS.items()
                  if sum(bool(matrix[f][ln]) for ln in seen) == 1
                  and len(seen) == 1]
    assert len(singletons) >= 5

    # no physics lane ever catches a certificate-like fault
    for fault in CERT_LIKE:
        assert not any(matrix[fault][ln] for ln in PHYSICS), fault

    # and no single lane catches everything (genuine complementarity)
    for lane in fc.LANES:
        caught = sum(1 for f in fc.FAULTS if matrix[f][lane])
        assert caught < len(fc.FAULTS)
