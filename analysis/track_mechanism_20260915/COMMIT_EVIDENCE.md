# EXEC_46 atomic commit line evidence

This ledger corrects the summary end-line typos in the bodies of `b215369`
and `b49276e` without rewriting history. The `git show --unified=0` hunks are
the authoritative before-to-after line evidence.

| Commit | File | Before -> after |
|---|---|---|
| `630ce62` | `include/EventAction.hh` | line 42 -> added line 43; line 60 -> added line 62 |
| `630ce62` | `src/RunAction.cc` | line 145 -> added line 146 |
| `630ce62` | `src/SiPMObservation.cc` | line 35 -> added line 37; line 42 -> added lines 45–54 |
| `630ce62` | `src/SiPMSD.cc` | line 89 -> added line 95; line 110 -> added line 117 |
| `19ebca6` | `src/RunAction.cc` | line 146 -> added lines 147–152 |
| `19ebca6` | `src/SiPMSD.cc` | line 78 -> replacement at lines 90–110; line 117 -> added lines 150–156 |
| `79bad4f` | `src/RunAction.cc` | line 152 -> added lines 153–154 |
| `79bad4f` | `src/SiPMSD.cc` | line 111 -> added lines 123–140; line 156 -> added lines 184–185 |
| `c8e40ce` | `include/PhotonTrackInfo.hh` | new file -> lines 1–17 |
| `c8e40ce` | `src/PhotonTrackInfo.cc` | new file -> lines 1–22 |
| `c8e40ce` | `src/RunAction.cc` | line 154 -> added line 155 |
| `c8e40ce` | `src/SteppingAction.cc` | line 160 -> added lines 163–180 |
| `c8e40ce` | `src/SiPMSD.cc` | line 139 -> added lines 142–153; line 185 -> added lines 200–201 |
| `4967ec8` | `src/RunAction.cc` | line 155 -> added line 156 |
| `4967ec8` | `src/SiPMSD.cc` | line 201 -> added lines 204–205 |
| `b215369` | `analysis/track_mechanism_20260915/analyze_validation.py` | new file -> lines 1–303 |
| `b49276e` | `analysis/track_mechanism_20260915/prepare_campaign.py` | new file -> lines 1–159 |

Reproduction command:

```bash
git show --format= --unified=0 <commit> --
```
