"""Smith-Waterman (local) and Needleman-Wunsch (global) sequence alignment."""


def _score(a: str, b: str, match: int, mismatch: int) -> int:
    return match if a == b else mismatch


def build_nw_matrix(
    seq1: str,
    seq2: str,
    match: int,
    mismatch: int,
    gap: int,
) -> list[list[int]]:
    """Build the Needleman-Wunsch (global) scoring matrix H."""
    rows, cols = len(seq1) + 1, len(seq2) + 1
    H = [[0] * cols for _ in range(rows)]

    for i in range(1, rows):
        H[i][0] = i * gap
    for j in range(1, cols):
        H[0][j] = j * gap

    for i in range(1, rows):
        for j in range(1, cols):
            diag = H[i - 1][j - 1] + _score(seq1[i - 1], seq2[j - 1], match, mismatch)
            up = H[i - 1][j] + gap
            left = H[i][j - 1] + gap
            H[i][j] = max(diag, up, left)
    return H


def build_sw_matrix(
    seq1: str,
    seq2: str,
    match: int,
    mismatch: int,
    gap: int,
) -> list[list[int]]:
    """Build the Smith-Waterman (local) scoring matrix H (zero-floored)."""
    rows, cols = len(seq1) + 1, len(seq2) + 1
    H = [[0] * cols for _ in range(rows)]

    for i in range(1, rows):
        for j in range(1, cols):
            diag = H[i - 1][j - 1] + _score(seq1[i - 1], seq2[j - 1], match, mismatch)
            up = H[i - 1][j] + gap
            left = H[i][j - 1] + gap
            H[i][j] = max(0, diag, up, left)
    return H


def global_align(
    seq1: str,
    seq2: str,
    match: int,
    mismatch: int,
    gap: int,
) -> tuple[str, str, int]:
    """Needleman-Wunsch alignment.

    Backtrace starts at the highest-scoring cell of the LAST COLUMN
    (semi-global on seq1's suffix) and ends at H[0][0].
    """
    H = build_nw_matrix(seq1, seq2, match, mismatch, gap)
    n = len(seq2)

    best_i = 0
    for i in range(len(H)):
        if H[i][n] > H[best_i][n]:
            best_i = i

    aligned1, aligned2 = [], []
    i, j = best_i, n
    score = H[best_i][n]

    while i > 0 or j > 0:
        current = H[i][j]
        if (
            i > 0
            and j > 0
            and current == H[i - 1][j - 1] + _score(seq1[i - 1], seq2[j - 1], match, mismatch)
        ):
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

    return "".join(reversed(aligned1)), "".join(reversed(aligned2)), score


def local_align(
    seq1: str,
    seq2: str,
    match: int,
    mismatch: int,
    gap: int,
) -> tuple[str, str, int]:
    """Smith-Waterman alignment.

    Backtrace starts at the highest-scoring cell anywhere in H and ends
    when H[i][j] == 0.
    """
    H = build_sw_matrix(seq1, seq2, match, mismatch, gap)

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
        diag = H[i - 1][j - 1]
        up = H[i - 1][j]
        left = H[i][j - 1]

        if current == diag + _score(seq1[i - 1], seq2[j - 1], match, mismatch):
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


def format_matrix(H: list[list[int]], seq1: str, seq2: str) -> str:
    """Render the NW matrix rotated: row 0 (gap-init, labeled U) at the bottom,
    seq1 chars labeling the rows above it, and a final header line `X U <seq2>`."""
    lines = []
    for i in range(len(H) - 1, -1, -1):
        label = "U" if i == 0 else seq1[i - 1]
        values = " ".join(str(v) for v in H[i])
        lines.append(f"{label} {values}")
    lines.append("X " + " ".join(["U"] + list(seq2)))
    return "\n".join(lines)
