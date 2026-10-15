# Magnetic-tweezers unzipping: data-processing pipeline (20261009 / 20261013)

Scripts, force calibration, DNA sequence and tether selections used for the dinucleosome force-ladder
(YTrunk_2601_980bp) unzipping data. Raw `.h5` / `.lut` data are **not** included.

## Method in one page
1. **Force**: magnet shift (mm) -> force with the `user_20261004` calibration
   (`scripts/process_tierA_sop.py: force_from_shift`; 0.70 / 1.00 / 1.70 mm = 46 / 30 / 18 pN, linear in log F,
   linear extrapolation outside). 5.00 mm ~ 1.6 pN (rest), 2.50 mm ~ 10 pN.
2. **Segmentation** (`process_20261009.py: segment_steps`): hold = shift within 0.004 mm of the step target; up / down moves in between.
3. **Smoothing** (`fec_all_20261009.py: smoothed_nm`): 5-frame running median, then 5-frame running mean; Z in nm.
4. **De-drift** (`fec_all_20261009.py: rest_drift`): median of the last 10 s of every 1.6 pN rest (>= 5 s long, >= 20 frames);
   points > 300 nm from the median of both neighbours dropped; piecewise-linear interpolation, held constant at the ends;
   fallback: line through the first and last second (< 2 good rests).
5. **Zero** (`fec_all_20261009.max_extension`): 3-MAD Hampel filter (+-1 s window), then the maximum of the 1 s running median of the 10 pN step.
6. **Survival** (`ladder_norm1028_20261009.bead_scale`): a step and all later ones are dropped if the bead is lost, flies (> 2000 nm),
   falls below -1000 nm, the next rest is not closed (> 300 nm from zero), or (last step, no rest after) rises > 600 nm in its last second.
7. **1028 reference**: same maximum on the highest surviving step; span = reference - zero; position = 1028 (ext - zero) / span.
8. **Absolute scale (model B)**: Luz_B = (ext - zero - 1800 bp handle stretch from 10 pN) / (2 x_ss(F));
   dsDNA 0.34 nm/bp, Lp 50 nm, S 1000 pN; ssDNA 0.58 nm/nt, Kuhn 1.55 nm, S 800 pN (`reprocess_26_core_sop.py`).
9. **Pooled peak maps** (`luzB_peakmap_by_force_20261009.py`): plateau samples from 0.3 s after the move, 5-unit bins, every tether weight 1,
   sqrt density, adaptive local-background peaks (basic-unzipping-analysis settings).

## Layout
- `scripts/` pipeline and analysis scripts (`process_*`, `fec_all_*`, `ladder_norm1028_*`, `state_occupancy_*`, `stall_*`, `smalljumps_*`, `*_20261013`).
  `whz_plot.py` + `assets/default.mplstyle` = figure helpers.
- `dna/` trunk sequence (971 bp amplicon, 5'->3' as in the design; the fork enters at the 3WJ end), design report, GenBank / FASTA, SnapGene map.
  601 sites in the design file: 304-450 and 501-647 (+4 bp stub offset: 314-460 / 511-657 in the analysis coordinates);
  in fork coordinates (30 bp T-arm included): first 601 u ~334-480 (dyad 407), second u ~531-677 (dyad 604).
- `selections/` the hand-selected tethers (data_id, ROI) of 20261009 (96) and 20261013 (101, +1 mM Mg2+).

## Running
Python 3.9+, numpy, pandas, scipy, h5py, matplotlib, PyMuPDF (fitz). The scripts use **absolute paths of the original machine**
(`/Users/whz/Desktop/DATA MT/unzipping/<date>` for raw data, `.../Proecessed Traces/WHZ_<date>` for output);
change `RAW` / `DEST` in `process_20261009.py` (batch scripts re-point them before importing) and the
`.venv` path before running. Batch 20261013 is run by re-pointing `process_20261009.RAW/DEST` (see `select_pages_20261013.py`).
