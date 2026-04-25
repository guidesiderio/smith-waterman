"""Tests for alignment.py."""

from alignment import build_nw_matrix, global_align, local_align


def test_nw_matrix_corner():
    """H[m][n] for ATC vs TCG with match=1, mismatch=-1, gap=-2 must be -2."""
    H = build_nw_matrix("ATC", "TCG", match=1, mismatch=-1, gap=-2)
    assert H[3][3] == -2, f"expected H[3][3] == -2, got {H[3][3]}"


def test_sw_score():
    """Best local score for ACACACTA vs AGCACACA must be 12 (match=2, mismatch=-1, gap=-1)."""
    _, _, score = local_align("ACACACTA", "AGCACACA", match=2, mismatch=-1, gap=-1)
    assert score == 12, f"expected score 12, got {score}"


def test_nw_backtrace_uses_last_column_max():
    """Global alignment score must equal the max of the last column, not H[m][n]."""
    seq1, seq2 = "ACGT", "AG"
    H = build_nw_matrix(seq1, seq2, match=2, mismatch=-1, gap=-1)
    n = len(seq2)
    last_col = [H[i][n] for i in range(len(H))]
    max_last = max(last_col)
    corner = H[len(seq1)][n]
    assert max_last != corner, "test setup invalid: max of last column equals H[m][n]"
    _, _, score = global_align(seq1, seq2, match=2, mismatch=-1, gap=-1)
    assert score == max_last, f"expected score {max_last}, got {score}"


def test_sw_alignment_non_empty():
    """Local alignment of related sequences must yield non-empty strings and a positive score."""
    a1, a2, score = local_align("GATTACA", "GCATGCU", match=2, mismatch=-1, gap=-1)
    assert score > 0, f"expected positive score, got {score}"
    assert a1 != "", "aligned seq1 should not be empty"
    assert a2 != "", "aligned seq2 should not be empty"


if __name__ == "__main__":
    test_nw_matrix_corner()
    test_sw_score()
    test_nw_backtrace_uses_last_column_max()
    test_sw_alignment_non_empty()
    print("All tests passed.")
