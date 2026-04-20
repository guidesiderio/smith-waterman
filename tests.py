"""Tests for the Smith-Waterman local sequence alignment implementation."""

from smith_waterman import smith_waterman


def test_classic_example():
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


def test_identical_sequences():
    """Identical sequences should align perfectly with no gaps.

    Every position matches, so score = match * len(seq).
    """
    seq = "GATTACA"
    a1, a2, score = smith_waterman(seq, seq, match=2, mismatch=-1, gap=-1)
    assert score == 2 * len(seq), f"expected score {2 * len(seq)}, got {score}"
    assert a1 == seq, f"expected '{seq}', got '{a1}'"
    assert a2 == seq, f"expected '{seq}', got '{a2}'"


def test_no_similarity():
    """Sequences with no matching characters should yield score 0 and empty alignment."""
    a1, a2, score = smith_waterman("AAAA", "CCCC", match=2, mismatch=-1, gap=-1)
    assert score == 0, f"expected score 0, got {score}"
    assert a1 == "", f"expected empty alignment, got '{a1}'"
    assert a2 == "", f"expected empty alignment, got '{a2}'"


if __name__ == "__main__":
    test_classic_example()
    test_identical_sequences()
    test_no_similarity()
    print("All tests passed.")
