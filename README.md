# Sequence Alignment Algorithms

Python implementations of two classic dynamic-programming algorithms for biological sequence alignment:

- **Smith-Waterman** — optimal *local* alignment
- **Needleman-Wunsch** — optimal *global* alignment

## Comparison

| Property | Smith-Waterman (local) | Needleman-Wunsch (global) |
|---|---|---|
| Alignment scope | Best-matching subsequence | Full sequences end-to-end |
| Matrix initialization | All zeros | First row/col filled with gap penalties |
| Cell minimum | 0 (no negatives) | Negatives allowed |
| Traceback start | Highest-scoring cell | Bottom-right corner H[m][n] |
| Traceback end | Cell reaches 0 | Top-left corner H[0][0] |
| Typical use | Finding conserved domains | Comparing closely related sequences |

---

## Smith-Waterman

### Usage

```python
from smith_waterman import smith_waterman

aligned1, aligned2, score = smith_waterman("ACACACTA", "AGCACACA")
print(score)     # 12
print(aligned1)  # A-CACACTA
print(aligned2)  # AGCACAC-A
```

### Parameters

| Parameter  | Default | Description                   |
|------------|---------|-------------------------------|
| `match`    | `2`     | Score for matching characters  |
| `mismatch` | `-1`    | Penalty for mismatches        |
| `gap`      | `-1`    | Penalty for gaps              |

### Return value

`(aligned_seq1, aligned_seq2, score)` — both aligned strings with `-` for gaps, and the optimal local alignment score.

---

## Needleman-Wunsch

### Usage

```python
from needleman_wunsch import needleman_wunsch

aligned1, aligned2, score = needleman_wunsch("ACACACTA", "AGCACACA")
print(score)     # 12
print(aligned1)  # A-CACACTA
print(aligned2)  # AGCACAC-A
```

### Parameters

| Parameter  | Default | Description                   |
|------------|---------|-------------------------------|
| `match`    | `2`     | Score for matching characters  |
| `mismatch` | `-1`    | Penalty for mismatches        |
| `gap`      | `-1`    | Penalty for gaps              |

### Return value

`(aligned_seq1, aligned_seq2, score)` — both aligned strings with `-` for gaps, and the optimal global alignment score.

---

## Running tests

```bash
python tests.py
```
