# Network-targeted TMS–EEG peak extraction

MATLAB code for semi-automatic extraction of early TMS-evoked potential
(TEP) peaks from subject-level grand-averaged TMS–EEG data, following a
collapsed-localizer strategy seeded on the grand-grand-average across all
participants and diagnostic groups.

This repository contains the peak-extraction script used in:

> Bertazzoli G., Canu E., Bagattini C., et al.
> *Network-targeted neurophysiological biomarkers of dysconnectivity and
> cognitive decline in Alzheimer's disease.*
> [https://doi.org/10.5281/zenodo.17162745](https://doi.org/10.5281/zenodo.17162745)


For the upstream preprocessing pipeline (SOUND + SSP-SIR + ICA), MRI and
DWI analyses, and downstream statistical analyses in R, please refer to
the corresponding sections of the paper and to the data repository linked
below.

---

## Procedure overview

Peak extraction follows the collapsed-localizer logic (Luck, 2014). The
operator's role is restricted to accepting, rejecting, or correcting an
algorithmically detected pick, within a window and ROI that are fixed
at the group level *before* subject-level extraction.

1. **Group-level definition of search window and ROI.** For each peak of
   interest, the time window and four-electrode ROI are defined on the
   grand-grand-average across all participants and diagnostic groups, and
   stored in `peak_extraction_grand_grand_average.xlsx`. Both are
   identical across subjects and independent of diagnostic group.

2. **Automatic candidate detection per subject.** For each subject and
   each peak, MATLAB's `findpeaks` is run on every ROI electrode within
   the predefined time window, with `MinPeakProminence = 0.05`. For
   negative peaks the signal is inverted prior to detection. The
   largest-prominence candidate across all ROI electrodes is selected as
   the automatic pick.

3. **Human accept / reject / override.** The TEP and the automatic pick
   are plotted in a maximized figure. The operator decides:
   - press **`1`** → accept the automatic pick;
   - press **`0`** → reject the peak (no peak present; latency is set
     to empty and the entry is flagged as missing downstream);
   - click on the trace, then press **`1`** → override the automatic
     pick with the manually clicked coordinates (used to correct
     clearly artifactual picks).

   The operator was blinded to diagnostic group at this step: figure
   titles and filenames contain only anonymized BIDS-style subject IDs,
   which do not encode diagnosis.

4. **Persistence and resume.** Results are written to disk after every
   accept/reject decision, so the review can be paused and resumed. Each
   inspection figure is saved as `.fig` and `.tif` for auditability.

5. **Long-to-wide reshape.** At the end of the script, the per-event
   struct (one entry per subject × peak) is reshaped into a wide
   subject-by-measure table and exported as CSV for downstream
   statistical analysis in R.

> [https://doi.org/10.5281/zenodo.17162745](https://doi.org/10.5281/zenodo.17162745)

---

## Outputs

All outputs are written to `./peak_extraction/`:

- `first_peak_pos_grand_gavg_avg_TP9_TP10_test.mat` — long-format struct
  array: one entry per (subject × peak) with fields `name`, `ampli`,
  `latency`, `channel`, `peakName`.
- `first_peak_pos_grand_gavg_avg_TP9_TP10_test.csv` — wide-format table:
  one row per subject; for each (peak × stimulation area), three columns
  (`amp`, `lat`, `elec`). Used as input for the R analysis scripts.
- One `.fig` + one `.tif` per (subject × peak), saving the trace and the
  accepted pick for visual auditability.



---

## Citation

If you use or adapt this code, please cite:

> Bertazzoli G., Canu E., Bagattini C., et al.
> *Network-targeted neurophysiological biomarkers of dysconnectivity and
> cognitive decline in Alzheimer's disease.*
> [https://doi.org/10.5281/zenodo.17162745](https://doi.org/10.5281/zenodo.17162745)

A BibTeX entry will be added here once the paper is published.

---

## Contact

For questions about the code or the underlying analyses:

- Giacomo Bertazzoli — `gbertazz@bidmc.harvard.edu`
- Marta Bortoletto — `marta.bortoletto@imtlucca.it`

---

## License

\[CC-BY-4.0]
