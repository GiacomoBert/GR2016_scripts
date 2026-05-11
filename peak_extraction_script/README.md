# Network-targeted TMS–EEG peak extraction

MATLAB code for semi-automatic extraction of early TMS-evoked potential (TEP) peaks from subject-level grand-averaged TMS–EEG data, using a collapsed-localizer strategy seeded on the grand-grand-average.

Companion code for:

> Bertazzoli G., Canu E., Bagattini C., et al. *Network-targeted neurophysiological biomarkers of dysconnectivity and cognitive decline in Alzheimer's disease.*
> [https://doi.org/10.5281/zenodo.17162745](https://doi.org/10.5281/zenodo.17162745)

---

## How it works

1. Time window and ROI electrodes are defined once on the grand-grand-average (all subjects, all groups) and stored in `peak_extraction_grand_grand_average.xlsx`. They are identical across subjects and independent of diagnostic group.
2. For each subject and peak, `findpeaks` runs on each ROI electrode within the predefined window (`MinPeakProminence = 0.05`; signal inverted for negative peaks). The largest-prominence candidate across electrodes is the automatic pick.
3. The trace and the pick are shown to the operator:
   - **`1`** → accept the automatic pick
   - **`0`** → reject (no peak; latency emptied)
   - click on trace, then **`1`** → override with the clicked coordinates
4. Figure titles and filenames carry only anonymized BIDS-style subject IDs, so the operator is blinded to diagnosis at this step.
5. Results are saved after every decision (resumable). At the end, the per-event struct is reshaped into a wide subject × measure CSV.

**Canonical search windows and ROI electrodes used in the published analysis are those reported in the paper (Table 2 / Methods) and on Zenodo: [https://doi.org/10.5281/zenodo.17162745](https://doi.org/10.5281/zenodo.17162745). Edit the local spreadsheet to match these before running.**

---

## Dependencies

- MATLAB R2020a+ (tested on R2024b)
- [EEGLAB](https://sccn.ucsd.edu/eeglab/) (`eeglab2024.2`)
- [FieldTrip](https://www.fieldtriptoolbox.org/) (`fieldtrip-20240731`)
- [`natsortfiles`](https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort)

Update the `addpath` calls and the absolute paths in the script to match your machine.

---

## Outputs

Written to `./peak_extraction/`:

- `first_peak_pos_grand_gavg_avg_TP9_TP10_test.mat` — long-format struct, one entry per subject × peak.
- `first_peak_pos_grand_gavg_avg_TP9_TP10_test.csv` — wide-format table, one row per subject. Input for the R analysis scripts.
- One `.fig` + one `.tif` per subject × peak.

---

## Citation

> Bertazzoli G., Canu E., Bagattini C., et al. *Network-targeted neurophysiological biomarkers of dysconnectivity and cognitive decline in Alzheimer's disease.*
> [https://doi.org/10.5281/zenodo.17162745](https://doi.org/10.5281/zenodo.17162745)

---

## Contact

- Giacomo Bertazzoli — `gbertazz@bidmc.harvard.edu`
- Marta Bortoletto — `marta.bortoletto@imtlucca.it`

---

## License

CC-BY-4.0
