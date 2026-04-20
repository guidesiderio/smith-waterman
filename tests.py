"""Tests for Smith-Waterman and Needleman-Wunsch sequence alignment implementations."""

from smith_waterman import smith_waterman
from needleman_wunsch import needleman_wunsch


def test_sw_classic_example():
    """Classic textbook example: ACACACTA vs AGCACACA.

    Expected optimal local alignment (match=2, mismatch=-1, gap=-1):
        A-CACACTA
        AGCACAC-A
    Score: 12
    """
    a1, a2, score = smith_waterman("ACACACTA", "AGCACACA")
    assert score == 12, f"expected score 12, got {score}"
    assert a1 == "A-CACACTA", f"expected 'A-CACACTA', got '{a1}'"
    assert a2 == "AGCACAC-A", f"expected 'AGCACAC-A', got '{a2}'"


def test_sw_identical_sequences():
    """Identical sequences should align perfectly with no gaps.

    Every position matches, so score = match * len(seq).
    """
    seq = "GATTACA"
    a1, a2, score = smith_waterman(seq, seq, match=2, mismatch=-1, gap=-1)
    assert score == 2 * len(seq), f"expected score {2 * len(seq)}, got {score}"
    assert a1 == seq, f"expected '{seq}', got '{a1}'"
    assert a2 == seq, f"expected '{seq}', got '{a2}'"


def test_sw_no_similarity():
    """Sequences with no matching characters should yield score 0 and empty alignment."""
    a1, a2, score = smith_waterman("AAAA", "CCCC", match=2, mismatch=-1, gap=-1)
    assert score == 0, f"expected score 0, got {score}"
    assert a1 == "", f"expected empty alignment, got '{a1}'"
    assert a2 == "", f"expected empty alignment, got '{a2}'"


def test_nw_classic_example():
    """Classic textbook example: ACACACTA vs AGCACACA.

    Expected optimal global alignment (match=2, mismatch=-1, gap=-1):
        A-CACACTA
        AGCACAC-A
    Score: 12
    """
    a1, a2, score = needleman_wunsch("ACACACTA", "AGCACACA")
    assert score == 12, f"expected score 12, got {score}"
    assert a1 == "A-CACACTA", f"expected 'A-CACACTA', got '{a1}'"
    assert a2 == "AGCACAC-A", f"expected 'AGCACAC-A', got '{a2}'"


def test_nw_identical_sequences():
    """Identical sequences should align perfectly with score = match * len(seq)."""
    seq = "GATTACA"
    a1, a2, score = needleman_wunsch(seq, seq, match=2, mismatch=-1, gap=-1)
    assert score == 2 * len(seq), f"expected score {2 * len(seq)}, got {score}"
    assert a1 == seq, f"expected '{seq}', got '{a1}'"
    assert a2 == seq, f"expected '{seq}', got '{a2}'"


def test_nw_one_gap_forced():
    """ACGT vs AGT forces one gap in the alignment.

    Expected global alignment (match=2, mismatch=-1, gap=-1):
        ACGT
        A-GT
    Score: 5
    """
    a1, a2, score = needleman_wunsch("ACGT", "AGT")
    assert score == 5, f"expected score 5, got {score}"
    assert a1 == "ACGT", f"expected 'ACGT', got '{a1}'"
    assert a2 == "A-GT", f"expected 'A-GT', got '{a2}'"


if __name__ == "__main__":
    test_sw_classic_example()
    test_sw_identical_sequences()
    test_sw_no_similarity()
    test_nw_classic_example()
    test_nw_identical_sequences()
    test_nw_one_gap_forced()
    print("All tests passed.")
