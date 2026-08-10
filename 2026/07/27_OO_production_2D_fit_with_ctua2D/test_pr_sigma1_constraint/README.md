# Prompt-sigma1 constrained lifetime-fit test

Scope: `1.0 < pT < 2.0 GeV`, `1.6 < |y| < 2.4`; all files here are isolated from the production directory.

The prompt test applies a Gaussian external constraint to `ctauTime1Scale`:

`sigma1 = 0.06859542 +/- 0.00823145 mm`.

`ctau_pr.C` and `ctau_np.C` print a pull summary after every fit: mean, RMS, counts above 2 and 3 standard deviations, and the maximum absolute pull.

Results are recorded in `pull_comparison.csv` and full logs/output snapshots are under `results/`.

Key finding: the sigma1 constraint is stable but does not change the nonprompt pull. The original nonprompt model had two log-lifetime parameters at their bounds. `ctau_np_relaxed2.C` relaxes those limits only in this test copy and improves chi2/ndf from 3.322 to 1.521.

Reproduce the final test from this directory:

```bash
root -l -b -q 'ctau_pr.C++(1.0,2.0,1.6,2.4,false,false,true,0.0685954216083519,0.00823145059300222)'
root -l -b -q 'ctau_np_relaxed2.C++(1.0,2.0,1.6,2.4,false,false)'
```
