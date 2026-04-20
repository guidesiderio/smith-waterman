# Smith-Waterman Local Sequence Alignment

Python implementation of the [Smith-Waterman algorithm](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) for optimal local alignment of two biological sequences.

## Usage

```python
from smith_waterman import smith_waterman

aligned1, aligned2, score = smith_waterman("ACACACTA", "AGCACACA")
print(score)     # 12
print(aligned1)  # A-CACACTA
print(aligned2)  # AGCACAC-A
```

### Parameters

| Parameter  | Default | Description                  |
|------------|---------|------------------------------|
| `match`    | `2`     | Score for matching characters |
| `mismatch` | `-1`    | Penalty for mismatches       |
| `gap`      | `-1`    | Penalty for gaps             |

### Return value

`(aligned_seq1, aligned_seq2, score)` — both aligned strings with `-` for gaps, and the optimal local alignment score.

## Running tests

```bash
python tests.py
```
