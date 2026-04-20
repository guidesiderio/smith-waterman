"""Smith-Waterman local sequence alignment algorithm."""

import itertools


def smith_waterman(
    seq1: str,
    seq2: str,
    match: int = 2,
    mismatch: int = -1,
    gap: int = -1,
) -> tuple[str, str, int]:
    """Align two sequences locally using the Smith-Waterman algorithm.

    Args:
        seq1: First sequence (rows).
        seq2: Second sequence (columns).
        match: Score awarded for matching characters.
        mismatch: Penalty for mismatching characters.
        gap: Penalty for a gap.

    Returns:
        Tuple of (aligned_seq1, aligned_seq2, score) where score is the
        optimal local alignment score.
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
    """Fill the Smith-Waterman scoring matrix H."""
    H = [[0] * cols for _ in range(rows)]

    for i in range(1, rows):
        for j in range(1, cols):
            diagonal = H[i - 1][j - 1] + _score(seq1[i - 1], seq2[j - 1], match, mismatch)
            up = H[i - 1][j] + gap
            left = H[i][j - 1] + gap
            H[i][j] = max(0, diagonal, up, left)

    return H


def _traceback(
    seq1: str,
    seq2: str,
    H: list[list[int]],
    match: int,
    mismatch: int,
    gap: int,
) -> tuple[str, str, int]:
    """Trace back from the highest-scoring cell to reconstruct the alignment."""
    best_score = 0
    best_i, best_j = 0, 0

    for i in range(len(H)):
        for j in range(len(H[0])):
            if H[i][j] > best_score:
                best_score = H[i][j]
                best_i, best_j = i, j

    aligned1, aligned2 = [], []
    i, j = best_i, best_j

    while i > 0 and j > 0 and H[i][j] > 0:
        current = H[i][j]
        diagonal = H[i - 1][j - 1]
        up = H[i - 1][j]
        left = H[i][j - 1]

        if current == diagonal + _score(seq1[i - 1], seq2[j - 1], match, mismatch):
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif current == up + gap:
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
        elif current == left + gap:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1
        else:
            raise ValueError(f"Traceback inconsistency at ({i}, {j})")

    return "".join(reversed(aligned1)), "".join(reversed(aligned2)), best_score


def _format_alignment(seq1: str, seq2: str) -> str:
    """Return a three-line alignment string with match indicators."""
    middle = "".join("|" if a == b else " " for a, b in itertools.zip_longest(seq1, seq2, fillvalue=" "))
    return f"{seq1}\n{middle}\n{seq2}"


if __name__ == "__main__":
    examples = [
        ("ACACACTA", "AGCACACA"),
        ("TGTTACGG", "GGTTGACTA"),
        ("AGTACGCA", "TATGC"),
    ]

    for s1, s2 in examples:
        a1, a2, score = smith_waterman(s1, s2)
        print(f"seq1: {s1}")
        print(f"seq2: {s2}")
        print(f"score: {score}")
        print(_format_alignment(a1, a2))
        print()
