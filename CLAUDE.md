# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

Run all tests:
```
python tests.py
```

Run with pytest (if installed):
```
pytest tests.py -v
```

Run the algorithm directly with built-in examples:
```
python smith_waterman.py
```

## Architecture

Single-module implementation in `smith_waterman.py`. The public API is one function:

```python
aligned1, aligned2, score = smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-1)
```

Internally split into:
- `_build_matrix` — fills the H scoring matrix (O(m×n) time and space)
- `_traceback` — finds the highest cell and walks back to reconstruct both aligned strings
- `_score` — match/mismatch helper

`tests.py` contains standalone test functions runnable as `python tests.py` (no test framework required). Note: the return order is `(aligned1, aligned2, score)` — not `(score, aligned1, aligned2)`.
