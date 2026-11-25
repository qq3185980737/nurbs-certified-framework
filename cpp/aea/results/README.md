# Results

- `e11_bench_summary.tsv`: one line per model from `data/run_bench.sh`; columns are
  `case-dir`, `case-file`, `RESULT ...`. The E11 table is the aggregate of lines
  containing `RESULT` for `fusion` and `mfcad`.
- `e10_abc_occt_pairs.txt`: face-level OCCT section results for the 43-face ABC model.
- `abc_bench.log`: full E10 benchmark log for the ABC model.

The current aggregate (τ = 0.02) is reproduced by:

```bash
python data/summarize.py --input results/e11_bench_summary.tsv
```
