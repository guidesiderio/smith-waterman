# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

Run the program (reads `input.txt`, writes `output.txt`):

```
python main.py
```

Run tests:

```
python tests.py
```

## Architecture

Two-file project plus an input file:

- `main.py` — entry point. Parses `input.txt`, calls into `alignment.py`,
  prints the formatted output and writes it to `output.txt`.
- `alignment.py` — core algorithms. Public functions:
  - `build_nw_matrix(seq1, seq2, match, mismatch, gap)` — Needleman-Wunsch matrix
  - `build_sw_matrix(seq1, seq2, match, mismatch, gap)` — Smith-Waterman matrix
  - `global_align(seq1, seq2, match, mismatch, gap)` — `(aligned1, aligned2, score)`
  - `local_align(seq1, seq2, match, mismatch, gap)` — `(aligned1, aligned2, score)`
  - `format_matrix(H, seq1, seq2)` — rotated matrix display
- `input.txt` — five lines: seq1, seq2, gap, mismatch, match.
- `tests.py` — single `__main__` block; assert-based.

## Backtrace conventions (assignment-specific)

- **Global (NW)** — backtrace starts at the cell with the highest score in the
  **last column** (semi-global on seq1's suffix), not at `H[m][n]`. Ends at `H[0][0]`.
- **Local (SW)** — backtrace starts at the highest-scoring cell anywhere in
  the matrix and stops when `H[i][j] == 0`.

## Output format

`main.py` prints (and writes to `output.txt`):

1. `# ** matrix **` header
2. The NW matrix, rendered rotated (row 0 labeled `U` at the bottom; header line `X U <seq2>` last)
3. A separator line of `=` characters
4. The score-parameter line: `** Match = X | mismatch = Y | Gap = Z **`
5. `Alinhamento Global` block: aligned seq1, aligned seq2 (chars space-separated), `Score = ...`
6. `Alinhamento Local` block: same shape

Return order from both align functions is `(aligned1, aligned2, score)`.
