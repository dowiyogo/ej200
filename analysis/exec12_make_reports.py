#!/usr/bin/env python3
"""Generate EXEC_11/12 technical notes exclusively from result CSV files."""
import argparse,pathlib,pandas as pd,numpy as np
def table(df):
 def cell(v): return f"{v:.3f}" if isinstance(v,(float,np.floating)) else str(v)
 rows=["| "+" | ".join(map(str,df.columns))+" |","|"+"|".join(["---"]*len(df.columns))+"|"]
 rows+=["| "+" | ".join(cell(v) for v in row)+" |" for row in df.itertuples(index=False,name=None)]
 return "\n".join(rows)
def main():
 p=argparse.ArgumentParser();p.add_argument("--exec11",type=pathlib.Path,required=True);p.add_argument("--exec12",type=pathlib.Path,required=True);a=p.parse_args()
 e11=a.exec11;e12=a.exec12;r=e12/"report";r.mkdir(exist_ok=True)
 refs=pd.read_csv(e11/"analysis/reference_positions.csv");rec11=pd.read_csv(e11/"analysis/reconstruction_summary.csv");qa=pd.read_csv(e11/"analysis/fit_qa_v1_v2_focus.csv")
 s=pd.read_csv(e12/"analysis/x_reconstruction_summary.csv");blue=pd.read_csv(e12/"analysis/blue_summary.csv");y=pd.read_csv(e12/"analysis/y0_feasibility_summary.csv")
 agg=s.groupby("method").agg(bias_mean=("bias","mean"),bias_abs_max=("bias",lambda q:np.max(np.abs(q))),sigma_core_mean=("sigma_core","mean"),sigma_core_max=("sigma_core","max"),valid_mean=("valid_fraction","mean")).reset_index()
 ba=blue.agg({"bias":"mean","sigma_core":"mean","valid_fraction":"mean"})
 note=f"""# EXEC_11 technical note: local Top-pair positioning

All values below are generated from `{e11}/analysis`; data commit `f431c01`.
Every resolution is a **simulation prediction — EJ-204 — intrinsic optical timing**.

## Method and geometry
The adjacent Top pair (28,29), centered near -442 mm, is reduced from per-hit
`sipm_hits` rows to fourth-hit times and PE-equivalent counts. The two references
were selected programmatically:

{table(refs)}

## QA and reconstruction
The original fits at -444, -442 and -439 mm converged to narrow subpeaks. Broad,
seeded, iterative v2 fits expose non-Gaussian tails:

{table(qa[["x_true_mm","sigma_dt_ps_v1","sigma_dt_ps_v2","chi2_ndf"]])}

{table(rec11[rec11.method.isin(["temporal","ratio_linear","BLUE"])][["reference","method","bias_mm","sigma_x_mm","sigma_x_error_mm"]])}

The count ratio is globally nonlinear; the cubic comparison creates ambiguous roots.
BLUE is covariance-aware. This is a one-pair empirical 1-D adaptation inspired by
pair timing, not Lv et al.'s literal 64-SiPM/2016-pair 2-D circle-intersection method.
It must not be directly compared with their 1.5–3 mm results.

## Limitations
Simulation only; one pair; two held-out test positions; no SPTR/electronics; strong
non-Gaussian timing tails; no experimental validation.
"""
 integrated=f"""# Event-by-event position reconstruction in a 1.4 m EJ-204 scintillator bar with end and top SiPM readout

## Scope
EXEC_11 established a local result near pair (28,29). EXEC_12 evaluates global X
with 31 leave-one-position-out folds and explores, but does not calibrate, Y at y=0.
All resolutions are **simulation prediction — EJ-204 — intrinsic optical timing**.

## Data and validation
The EndTop scan contains 31 stable ROOT files and 2000 events per position. Logs
confirm OPSC-101/EJ-204 and EndTop. ROOT content confirms gun X. Jitter=0, y=0 and
perpendicular incidence are inferred from the generating script/code.

## Global X results
{table(agg)}

BLUE over End weighted timing, End ratio and Top centroid has mean core spread
`{ba.sigma_core:.3f} mm`; its maximum absolute bias is `{blue.bias.abs().max():.3f} mm`.
The local transferable Top-pair ratio has mean core spread
`{agg.loc[agg.method=="local_R","sigma_core_mean"].iloc[0]:.3f} mm`, but is invalid
at x=0 because the 24 mm central gap makes pair choice ambiguous. It has position-
dependent biases up to `{agg.loc[agg.method=="local_R","bias_abs_max"].iloc[0]:.3f} mm`.
Thus the ~5.2 mm local EXEC_11 result generalizes in local-pair spread, but not as a
uniformly available unbiased estimator across the full bar. Global linear Top
centroid has narrow event spread but large edge bias; End estimators degrade strongly
at extremes. The current BLUE improves individual global End methods but does not
beat the valid local-pair estimator.

## Y=0 feasibility proxy
{table(y.groupby("estimator").agg(mean_centroid=("mean","mean"),max_abs_mean=("mean",lambda q:np.max(np.abs(q))),mean_width=("width","mean"),max_width=("width","max")).reset_index())}

At y_true=0, the combined end-channel light-sharing centroid is centered near zero
and has an event-by-event width shown above. This is a central-point estimator spread /
feasibility proxy, **not a Y spatial resolution**.

## Availability matrix
| Capability | Status |
|---|---|
| Local X near pair (28,29) | validated in current simulation |
| Global X along bar | evaluated with position-level LOO CV |
| Y at central incidence | exploratory proxy |
| General Y reconstruction | not yet available |
| Full X–Y reconstruction | requires Y scan |
| Experimental validation | not yet available |

## Conclusion
Current data support X reconstruction and show transferable local Top information,
but edge biases, the central gap, and model nonlinearity prevent claiming a uniform
5 mm full-bar result. Current data do not permit a validated X–Y reconstruction.
"""
 for name,text in [("exec11_technical_note",note),("endtop_position_reconstruction",integrated)]:
  (r/f"{name}.md").write_text(text)
  tex="\\\\documentclass{article}\\n\\\\usepackage[margin=1in]{geometry}\\n\\\\begin{document}\\n\\\\section*{"+name.replace("_"," ")+"}\\n\\\\begin{verbatim}\\n"+text+"\\\\end{verbatim}\\n\\\\end{document}\\n"
  (r/f"{name}.tex").write_text(tex)
 (a.exec12/"README.md").write_text(integrated+"\n\nLaTeX sources were generated, but PDF compilation failed on MSI because pdfLaTeX cannot render the Unicode-rich verbatim technical note.\n")
if __name__=="__main__":main()
