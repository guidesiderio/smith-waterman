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

Run algorithms directly with built-in examples:
```
python smith_waterman.py
python needleman_wunsch.py
```

## Architecture

Two single-module implementations sharing the same public API signature:

```python
aligned1, aligned2, score = smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-1)
aligned1, aligned2, score = needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap=-1)
```

Each module is internally split into:
- `_build_matrix` — fills the H scoring matrix (O(m×n) time and space)
- `_traceback` — walks back through H to reconstruct both aligned strings
- `_score` — match/mismatch helper
- `_format_alignment` — three-line display string with match indicators

Key difference: Smith-Waterman (local) floors cells at 0 and traces back from the highest-scoring cell; Needleman-Wunsch (global) allows negative values and always traces back from H[m][n] to H[0][0].

`tests.py` contains standalone test functions runnable as `python tests.py` (no test framework required). Note: the return order is `(aligned1, aligned2, score)` — not `(score, aligned1, aligned2)`.
