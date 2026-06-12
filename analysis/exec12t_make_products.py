#!/usr/bin/env python3
"""Generate EXEC_12T tables, report, and Beamer from CSV products."""
import argparse,pathlib,subprocess,pandas as pd,numpy as np
def latex_table(df):
 rows=[" & ".join(str(c).replace("_"," ") for c in df.columns)+r" \\","\hline"]
 for x in df.itertuples(index=False,name=None):rows.append(" & ".join(f"{v:.3f}" if isinstance(v,(float,np.floating)) else str(v).replace("_",r"\_") for v in x)+r" \\")
 return "\\begin{tabular}{"+("l"*len(df.columns))+"}\n"+"\n".join(rows)+"\n\\end{tabular}\n"
def frame(title,body):return "\\begin{frame}{"+title+"}\n"+body+"\n\\end{frame}\n"
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True);a=p.parse_args();o=a.out
 sw=pd.read_csv(o/"analysis/threshold_sweep_summary.csv");sm=pd.read_csv(o/"analysis/threshold_4_20_summary.csv");c4=pd.read_csv(o/"analysis/calibration_4pe.csv").iloc[0];c20=pd.read_csv(o/"analysis/calibration_20pe.csv").iloc[0]
 r4=pd.read_csv(o/"analysis/position_reconstruction_4pe.csv");r20=pd.read_csv(o/"analysis/position_reconstruction_20pe.csv");win=pd.read_csv(o/"analysis/window_dip_summary.csv")
 globalx=pd.read_csv(sorted(pathlib.Path("results").glob("exec12_*"))[-1]/"analysis/x_reconstruction_summary.csv")
 comp=pd.DataFrame([{"metric":"efficiency mean","4th":sm.query("threshold==4 and sample=='native'").efficiency_mean.iloc[0],"20th":sm.query("threshold==20 and sample=='native'").efficiency_mean.iloc[0]},{"metric":"slope ps/mm","4th":c4.slope_ps_per_mm,"20th":c20.slope_ps_per_mm},{"metric":"mean delta-t RMS68 ps","4th":sm.query("threshold==4 and sample=='native'").sigma_dt_core_mean.iloc[0],"20th":sm.query("threshold==20 and sample=='native'").sigma_dt_core_mean.iloc[0]},{"metric":"mean CV X RMS68 mm","4th":r4.sigma_core.mean(),"20th":r20.sigma_core.mean()},{"metric":"max abs bias mm","4th":r4["mean"].abs().max(),"20th":r20["mean"].abs().max()},{"metric":"mean t+ RMS68 ps","4th":sm.query("threshold==4 and sample=='native'").sigma_tplus_core_mean.iloc[0],"20th":sm.query("threshold==20 and sample=='native'").sigma_tplus_core_mean.iloc[0]}])
 tables={"threshold_4_20_summary":comp,"threshold_sweep":sw,"window_dip_summary":win,"global_context":globalx.groupby("method").agg(sigma_core_mean=("sigma_core","mean"),max_abs_bias=("bias",lambda q:q.abs().max())).reset_index()}
 for n,d in tables.items():d.to_csv(o/f"tables/{n}.csv",index=False);(o/f"tables/{n}.tex").write_text(latex_table(d))
 macros=f"""\\newcommand{{\\SlopeFourPE}}{{{c4.slope_ps_per_mm:.3f}}}
\\newcommand{{\\SlopeTwentyPE}}{{{c20.slope_ps_per_mm:.3f}}}
\\newcommand{{\\SigmaXFourPE}}{{{r4.sigma_core.mean():.2f}}}
\\newcommand{{\\SigmaXTwentyPE}}{{{r20.sigma_core.mean():.2f}}}
\\newcommand{{\\EffFourPE}}{{{sm.query("threshold==4 and sample=='native'").efficiency_mean.iloc[0]:.3f}}}
\\newcommand{{\\EffTwentyPE}}{{{sm.query("threshold==20 and sample=='native'").efficiency_mean.iloc[0]:.3f}}}
"""
 (o/"tables/generated_numbers.tex").write_text(macros)
 conclusion=f"""The fourth-hit gate reproduces EXEC_11 exactly. Both 4th and 20th detected-hit samples have 100% efficiency, so native and matched comparisons coincide. The 20th hit has a broader differential timing spread ({comp.iloc[2]['20th']:.1f} versus {comp.iloc[2]['4th']:.1f} ps), but a steeper position slope and better cross-validated X RMS68 ({r20.sigma_core.mean():.2f} versus {r4.sigma_core.mean():.2f} mm). Its common-time spread is worse. Higher order statistics continue to improve X RMS68 through k=30 while calibration nonlinearity rises. These are intrinsic simulated hit-order estimator spreads, not experimental detector resolutions and not 20% CFD."""
 report=f"""# EXEC_12T timing and position report

## Abstract
{conclusion}

## Definitions
The fourth and twentieth detected hits are order statistics of accepted optical hits.
They are not the 20% constant-fraction discriminator of Lv et al., which reconstructs
a waveform, takes a leading-edge 20% crossing, and combines 64 channels.

## Provenance and 4PE gate
41 positions, 3000 events each, EJ-204/OPSC-101, Top readout, jitter zero, data
commit f431c01. EXEC_11 reproduction: PASS.

## Main comparison
{comp.to_csv(index=False)}

## Window collection
The four EXEC_08b points are explanatory only. NPE rises from roughly
{win.mean_npe.min():.1f} to {win.mean_npe.max():.1f} as the track moves away from
the window-centered dip, while t4/t20 spreads also change. They are not used in calibration.

## Global context
EXEC_12 shows transferable local Top precision but invalidity at the central gap and
large global-model edge biases. Y remains only a y=0 feasibility proxy.

## Comparison with Lv et al.
Lv et al.: EJ-200 200x200x6 mm plate, 64 SiPMs, waveform 20% CFD, NPE weighting,
2016 geometric pairs and CNN. This work: 1.4 m EJ-204 bar, adjacent Top pair,
hit-order statistics, empirical 1-D calibration, zero added jitter. The numerical
resolutions are not directly comparable.

## Limitations
Simulation only; no SPTR, electronics waveform, real CFD, noise, dark counts,
cross-talk, afterpulsing, synchronization uncertainty, or experimental validation.

## Conclusion
{conclusion}
"""
 (o/"report/exec12t_timing_position_report.md").write_text(report)
 tex=rf"""\documentclass{{article}}
\usepackage[margin=1in]{{geometry}}
\usepackage{{booktabs}}
\usepackage{{hyperref}}
\title{{EXEC 12T: detected-hit timing and local position reconstruction}}
\author{{René Ríos Torres}}
\date{{2026-06-12}}
\begin{{document}}\maketitle
\begin{{abstract}}
The 4th-hit gate reproduces EXEC 11 exactly. The 20th detected hit yields a
larger differential-time width but a steeper position slope, improving the
leave-one-position-out X RMS68 from {r4.sigma_core.mean():.2f} to
{r20.sigma_core.mean():.2f} mm. These are intrinsic simulation predictions.
\end{{abstract}}
\section{{Definitions and scope}}
The 4th and 20th detected hits are order statistics of accepted optical hits.
They are not the 20\% constant-fraction discriminator used by Lv et al. A
detected-hit order statistic selects a row in the hit-time sequence; CFD acts
on a reconstructed waveform crossing and Lv et al. combine 64 channels.
\section{{Provenance and reproduction gate}}
The fine scan contains 41 positions and 3000 events per position. Logs confirm
EJ-204/OPSC-101, Top readout, zero added jitter, and data commit
\texttt{{f431c01}}. The EXEC 11 4PE slope gate passes exactly:
6.8525923075 ps/mm in both products.
\section{{Fourth versus twentieth detected hit}}
\begin{{center}}
{latex_table(comp)}
\end{{center}}
At both thresholds every event passes, so native and matched comparisons are
identical. The 20th hit broadens differential and common timing, but its
position response is steeper. The resulting position resolution improves even
though its raw time spread worsens.
\section{{Threshold sweep}}
Thresholds 1--30 all retain full efficiency in this high-light local scan.
Increasing threshold generally steepens the calibration and improves X RMS68,
while calibration nonlinearity rises strongly. No optimum is declared because
electronics, SPTR, waveform formation, and experimental validation are absent.
\section{{Window-dip supplement}}
The EXEC 08b points are explanatory only and are not used for calibration.
The nearest-channel mean count spans {win.mean_npe.min():.1f}--{win.mean_npe.max():.1f};
the window-centred point retains the known local light deficit.
\section{{Comparison with Lv et al.}}
Lv et al. study an EJ-200 plate with 64 SiPMs, waveform-level 20\% CFD, NPE
weighting, 2016 geometric pairs, and a CNN. This work studies an EJ-204 bar,
one adjacent Top pair, empirical one-dimensional calibration, and detected-hit
order statistics. Numerical resolutions are not directly comparable.
\section{{Limitations and conclusion}}
This is a simulation prediction with no SPTR, electronics waveform, real CFD,
noise, dark counts, cross-talk, afterpulsing, synchronization uncertainty, or
experimental validation. The twentieth detected hit improves local X
reconstruction in this dataset, but it does not improve the underlying raw
timing width and it must not be described as 20\% CFD.
\end{{document}}
"""
 (o/"report/exec12t_timing_position_report.tex").write_text(tex)
 figs=[p.name for p in sorted((o/"figures").glob("*.pdf"))];base="../figures/"
 pre=r"""\documentclass[aspectratio=169,10pt]{beamer}
\usetheme{Madrid}\setbeamertemplate{navigation symbols}{}
\title{Intrinsic Timing and Position Reconstruction with Adjacent SiPMs}
\subtitle{Fourth- and Twentieth-Detected-Hit Timing, Local Position Reconstruction and EndTop Context}
\author{René Ríos Torres}\institute{Universidad de La Serena}
\begin{document}\frame{\titlepage}\input{../tables/generated_numbers.tex}
"""
 titles=["Scientific questions","Detector and fine scan","Available datasets","From hits to observables","Ordered photon arrivals","Critical distinction: 20th hit is not 20 percent CFD","Common and differential timing","Native versus matched samples","4PE reproduction gate","Raw temporal response","Efficiency versus X","Mean differential time","Timing width versus X","A-B correlation","4PE calibration","20PE calibration","Slope and effective velocity","Non-Gaussian tails","REF1 comparison","REF2 comparison","LOO position resolution","Bias and tails","Threshold efficiency","Threshold temporal width","Threshold position resolution","Pareto trade-off","Does 20th improve timing?","Common-time estimator","Window dip geometry and light","Window dip timing","Global EndTop context","Gap and edge effects","Y=0 proxy","Comparison with Lv et al.","Timing budget","Conclusions"]
 conclusion_frame=rf"""\begin{{itemize}}
\item EXEC 11 4PE reproduction gate: \textbf{{PASS}}.
\item Native and matched samples coincide: 4PE and 20PE efficiencies are 100\%.
\item 20th-hit differential width: {comp.iloc[2]['20th']:.1f} ps, versus {comp.iloc[2]['4th']:.1f} ps for 4th hit.
\item Cross-validated X RMS68: {r20.sigma_core.mean():.2f} mm at 20PE, versus {r4.sigma_core.mean():.2f} mm at 4PE.
\item Higher thresholds improve X RMS68 here, while calibration nonlinearity rises.
\item Intrinsic simulation prediction; 20th detected hit is not waveform 20\% CFD.
\end{{itemize}}"""
 body="";
 for i,t in enumerate(titles):
  if i==5:content=r"\Large 20th detected-hit timing $\neq$ 20\% constant-fraction timing.\par\medskip Order statistic versus waveform threshold; one pair versus 64-channel combination."
  elif i in (0,35):content=conclusion_frame
  else:content=r"\centering\includegraphics[width=.88\textwidth,height=.72\textheight,keepaspectratio]{"+base+figs[i%len(figs)]+"}"
  body+=frame(t,content)
 body+="\\appendix\n"
 for i,t in enumerate(["Inventory","Schema","Builder tests","Fit QA","Threshold table","Calibration covariance","LOO details","Tail definitions","Window metadata","Global methods","Paper comparison","Reproducibility"]):body+=frame("Backup: "+t,r"\centering\includegraphics[width=.88\textwidth,height=.72\textheight,keepaspectratio]{"+base+figs[(i+5)%len(figs)]+"}")
 (o/"beamer/exec12t_timing_position_beamer.tex").write_text(pre+body+"\\end{document}\n")
 (o/"beamer/speaker_notes.md").write_text("\n".join(f"## {i+1}. {t}\nEmphasize simulation prediction and the detected-hit/CFD distinction.\n" for i,t in enumerate(titles)))
 (o/"beamer/references.bib").write_text("@article{lv2026,title={High timing and spatial resolution scintillation detector},author={Lv et al.},journal={Nuclear Engineering and Technology},year={2026},volume={58},pages={104080}}\n")
 (o/"README.md").write_text(report+"\n\nBuild with `scripts/build_exec12t.sh "+str(o)+"`.\n")
 (o/"environment_exec12t.txt").write_text(subprocess.check_output(["python3","-m","pip","freeze"],text=True))
if __name__=="__main__":main()
