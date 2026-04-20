"""Needleman-Wunsch global sequence alignment algorithm."""

from itertools import zip_longest


def needleman_wunsch(
    seq1: str,
    seq2: str,
    match: int = 2,
    mismatch: int = -1,
    gap: int = -1,
) -> tuple[str, str, int]:
    """Align two sequences globally using the Needleman-Wunsch algorithm.

    Args:
        seq1: First sequence (rows).
        seq2: Second sequence (columns).
        match: Score awarded for matching characters.
        mismatch: Penalty for mismatching characters.
        gap: Penalty for a gap.

    Returns:
        Tuple of (aligned_seq1, aligned_seq2, score) where score is the
        optimal global alignment score.
    """
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    H = _build_matrix(seq1, seq2, rows, cols, match, mismatch, gap)
    aligned1, aligned2, score = _traceback(seq1, seq2, H, match, mismatch, gap)

    return aligned1, aligned2, score


def _score(a: str, b: str, match: int, mismatch: int) -> int:
    return match if a == b else mismatch


def _build_matrix(
    seq1: str,
    seq2: str,
    rows: int,
    cols: int,
    match: int,
    mismatch: int,
    gap: int,
) -> list[list[int]]:
    """Fill the Needleman-Wunsch scoring matrix H."""
    H = [[0] * cols for _ in range(rows)]

    for i in range(1, rows):
        H[i][0] = i * gap
    for j in range(1, cols):
        H[0][j] = j * gap

    for i in range(1, rows):
        for j in range(1, cols):
            diagonal = H[i - 1][j - 1] + _score(seq1[i - 1], seq2[j - 1], match, mismatch)
            up = H[i - 1][j] + gap
            left = H[i][j - 1] + gap
            H[i][j] = max(diagonal, up, left)

    return H


def _traceback(
    seq1: str,
    seq2: str,
    H: list[list[int]],
    match: int,
    mismatch: int,
    gap: int,
) -> tuple[str, str, int]:
    """Trace back from H[m][n] to H[0][0] to reconstruct the alignment."""
    i, j = len(seq1), len(seq2)
    best_score = H[i][j]

    aligned1, aligned2 = [], []

    while i > 0 or j > 0:
        current = H[i][j]

        if i > 0 and j > 0 and current == H[i - 1][j - 1] + _score(seq1[i - 1], seq2[j - 1], match, mismatch):
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and current == H[i - 1][j] + gap:
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1

    return "".join(reversed(aligned1)), "".join(reversed(aligned2)), best_score


def _format_alignment(seq1: str, seq2: str) -> str:
    """Return a three-line alignment string with match indicators."""
    middle = "".join("|" if a == b else " " for a, b in zip_longest(seq1, seq2, fillvalue=" "))
    return f"{seq1}\n{middle}\n{seq2}"


if __name__ == "__main__":
    examples = [
        ("ACACACTA", "AGCACACA"),
        ("TGTTACGG", "GGTTGACTA"),
        ("AGTACGCA", "TATGC"),
    ]

    for s1, s2 in examples:
        a1, a2, score = needleman_wunsch(s1, s2)
        print(f"seq1: {s1}")
        print(f"seq2: {s2}")
        print(f"score: {score}")
        print(_format_alignment(a1, a2))
        print()
