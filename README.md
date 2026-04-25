# Sequence Alignment (Needleman-Wunsch + Smith-Waterman)

Single-entry-point Python implementation of two classic dynamic-programming
sequence alignment algorithms:

- **Needleman-Wunsch** — global alignment
- **Smith-Waterman** — local alignment

The program reads parameters from `input.txt`, runs both alignments, prints the
result to the screen, and saves it to `output.txt`.

## Files

| File           | Purpose                                                        |
|----------------|----------------------------------------------------------------|
| `main.py`      | Entry point — reads `input.txt`, runs both, prints + saves     |
| `alignment.py` | Core logic — matrix building and traceback for NW and SW       |
| `input.txt`    | Input parameters (sequences + scores)                          |
| `output.txt`   | Generated when `main.py` runs                                  |
| `tests.py`     | Test suite                                                     |

## Input format (`input.txt`)

```
seq1            (vertical, rows)
seq2            (horizontal, columns)
gap penalty     (integer)
mismatch        (integer)
match           (integer)
```

Example:

```
ATC
TCG
-2
-1
1
```

## Running

```bash
python main.py
```

Example output for the input above:

```
# ** matrix **

C -6 -3 0 -2
T -4 -1 -2 -4
A -2 -1 -3 -5
U 0 -2 -4 -6
X U T C G
===========================================================
** Match = 1 | mismatch = -1 | Gap = -2 **

Alinhamento Global
A T C -
- T C G
Score = -2

Alinhamento Local
T C
T C
Score = 2
```

## Tests

```bash
python tests.py
```

## Algorithm details

The matrix shown is the **Needleman-Wunsch** scoring matrix. It is rendered
rotated: row 0 (gap-initialized, labeled `U`) appears at the bottom, with
`seq1` characters labeling the rows above. The header line `X U <seq2>` sits
beneath the data.

Backtrace rules:

- **Global (NW)** — full NW matrix; backtrace begins at the cell with the
  highest score in the **last column** (not the bottom-right corner) and
  ends at `H[0][0]`.
- **Local (SW)** — zero-floored matrix; backtrace begins at the highest
  cell anywhere in the matrix and stops when `H[i][j] == 0`.
