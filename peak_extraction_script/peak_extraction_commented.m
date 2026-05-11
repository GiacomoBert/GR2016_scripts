% =========================================================================
% peak_extraction.m
% -------------------------------------------------------------------------
% Semi-automatic extraction of early TMS-evoked potential (TEP) peaks from
% subject-level grand-averaged datasets, following a collapsed-localizer
% strategy seeded on the grand-grand-average.
%
% Companion code for:
%   Bertazzoli G., Canu E., Bagattini C., et al.
%   "Network-targeted neurophysiological biomarkers of dysconnectivity and
%   cognitive decline in Alzheimer's disease."
%   [journal, year, DOI when available]
%
% Author : Giacomo Bertazzoli (gbertazz@bidmc.harvard.edu)
% License: <add license here, e.g. MIT or CC-BY-4.0>
%
% -------------------------------------------------------------------------
% PROCEDURE OVERVIEW
% -------------------------------------------------------------------------
% Peak extraction follows the collapsed-localizer logic (Luck, 2014):
%
%   1. The TIME WINDOW and ELECTRODE ROI for each peak of interest are
%      defined once, on the grand-grand-average across all participants
%      and diagnostic groups. These are stored in an external Excel file
%      (see `peak_extraction_grand_grand_average.xlsx`). Both are fixed for
%      all subjects and are therefore independent of group membership.
%
%   2. For each subject and each peak, the algorithm runs MATLAB's
%      `findpeaks` on every electrode in the ROI, within the predefined
%      time window, with `MinPeakProminence = 0.05`. For negative peaks,
%      the signal is inverted prior to running `findpeaks`. The
%      largest-prominence candidate across all ROI electrodes is selected
%      as the automatic pick.
%
%   3. The TEP and the automatic pick are plotted in a maximized figure
%      and a human operator decides whether the pick is acceptable:
%
%        - press '1'  -> accept the automatic pick;
%        - press '0'  -> reject the peak (no peak considered present;
%                        latency is set to empty and the entry is
%                        flagged as missing downstream);
%        - click on the trace, then press '1' -> override the automatic
%                        pick with the manually clicked coordinates
%                        (used to correct clearly artifactual picks).
%
%      The operator was blinded to diagnostic group at this step: figure
%      titles and filenames contain only the anonymized BIDS-style subject
%      ID, which does not encode diagnosis.
%
%   4. Amplitude, latency, channel, and peak name are stored in a struct
%      and saved to disk after every accept/reject decision, so the
%      review can be paused and resumed (resume logic at the top of the
%      main loop, controlled by `overwright`).
%
%   5. At the end of the script, the per-event struct is reshaped into a
%      wide subject-by-measure table and exported as CSV for downstream
%      statistical analysis in R.
%
% -------------------------------------------------------------------------
% DEPENDENCIES
% -------------------------------------------------------------------------
%   - MATLAB R2020a or later (tested with R2024b)
%   - EEGLAB (tested with eeglab2024.2)
%   - FieldTrip (tested with fieldtrip-20240731)
%   - natsortfiles (File Exchange, Stephen Cobeldick)
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
%   - `data_all_avg.mat`
%       Cell array of subject-level FieldTrip timelock structures, one per
%       stimulation site (LF, RF, LP, RP), pre-averaged across trials.
%
%   - `peak_extraction_grand_grand_average.xlsx`, sheet 'TEP_test'
%       Table specifying, for each peak of interest:
%           area      : stimulation site (LF / RF / LP / RP)
%           peakName  : peak identifier (e.g. 'N20F_amp_L_DLPFC')
%           start     : start of the search window (seconds)
%           xEnd      : end of the search window (seconds)
%           elec      : space-separated list of ROI electrodes
%
%   - `TEP_DTI_measures_total_*.xlsx`, sheet 'Sheet1'
%       Master subject list (column `ID_subj_NT`) used to restrict the
%       analysis to participants included in the present study.
%
% -------------------------------------------------------------------------
% OUTPUTS (written to ./peak_extraction/)
% -------------------------------------------------------------------------
%   - first_peak_pos_grand_gavg_avg_TP9_TP10_test.mat
%       Long-format struct array: one entry per (subject x peak)
%       containing fields: name, ampli, latency, channel, peakName.
%
%   - first_peak_pos_grand_gavg_avg_TP9_TP10_test.csv
%       Wide-format table: one row per subject; one set of columns
%       (amplitude, latency, electrode) per (peak x stimulation area).
%       Imported by the R analysis scripts.
%
%   - One .fig and one .tif file per (subject x peak), saving the trace
%       and the accepted pick for visual auditability.
%
% -------------------------------------------------------------------------
% NOTES ON REPRODUCIBILITY
% -------------------------------------------------------------------------
%   - The search window and ROI are fixed at the group level prior to
%     subject-level extraction. The operator does not choose where or when
%     to look; the operator only accepts, rejects, or corrects the
%     algorithm's pick within the predefined window and ROI.
%
%   - All absolute paths below reflect the original analysis machine.
%     Adapt them to your local layout before running.
% =========================================================================


%% Reset workspace and define base directory
clear all
close all
clc

% Resolve the directory in which this script lives. Used as the working
% root for the output folder created below.
main_dir = extractBefore(mfilename( 'fullpath' ), ['\' mfilename]);


%% Add required plug-ins
% Paths below reflect the original analysis machine. Replace with the
% local paths to natsortfiles, EEGLAB and FieldTrip on your system.
addpath('E:\BrightFocus_PreclinicalAD\01_preprocessing\plugins\natsortfiles');   % natsort (File Exchange)
addpath('E:\BrightFocus_PreclinicalAD\01_preprocessing\plugins\eeglab2024.2');   % EEGLAB
addpath('E:\BrightFocus_PreclinicalAD\01_preprocessing\plugins\fieldtrip-20240731'); % FieldTrip


%% Create output folder
% All outputs (struct, CSV, per-peak figures) are written here.
savepath = [main_dir '/peak_extraction'];
mkdir(savepath);


%% Initialize EEGLAB and FieldTrip without launching the EEGLAB GUI
eeglab nogui
ft_defaults


%% Paths to raw and cleaned BIDS data (kept for reference; not used here)
% These pointers describe where the underlying data come from; they are
% recorded here for traceability but are not read by this script.
raw_folder  = 'F:\gbertazz\BERTANAS\archive\raw_data';
data_parent = 'F:\gbertazz\BERTANAS\archive\cleaned_data\derivatives';


%% Load subject-level grand averages
% `data_all_avg` is a 1xN cell array, one cell per stimulation site
% (LF, RF, LP, RP). Each cell contains a list of FieldTrip timelock
% structures, one per subject and session, already averaged across trials.
load('F:\gbertazz\BERTANAS\giacomo\disco_F\GR2016_20230119\04_analisi_stat\grand_averages\gavg\grand_gavg_avg_2024\data_all_avg.mat');


%% Load the master subject list and filter datasets
% The master spreadsheet `TEP_DTI_measures_total_*.xlsx` contains the IDs
% of every participant retained for the present study. Datasets in
% `data_all_avg` that do not match this list are removed before peak
% extraction.
filename   = 'H:\GR2016_peer\TEP_DTI_measures_total_20240315.xlsx';
sheet      = 'Sheet1';
data       = readtable(filename, 'sheet', sheet);
columnName = 'ID_subj_NT';                       % column holding subject IDs
subj       = table2cell(data(:, columnName));

% Strip diagnostic-group letter (if any) from each ID and keep only the
% leading 6-character anonymous identifier. This ensures filtering is by
% subject ID alone, with no group information embedded.
subj          = extractBefore(subj, 7);
filterStrings = subj;

% Remove any dataset whose name does not match the filtered subject list.
% For each stimulation site, empty out non-matching datasets and then
% drop the now-empty columns.
for x = 1:length(data_all_avg)
    for y = 1:length(data_all_avg{1,x})
        current_dataset_name = data_all_avg{1,x}{y}.setname;
        if ~contains(current_dataset_name, filterStrings)
            data_all_avg{1,x}{y} = [];
        end
    end
    % Detect empty cells
    isEmptyElement    = cellfun(@isempty, data_all_avg{1,x});
    % Count empties per column
    emptyColumnsCount = sum(isEmptyElement, 1);
    % Keep only columns that are not entirely empty
    columnsToKeep     = emptyColumnsCount ~= size(data_all_avg{1,x}, 1);
    data_all_avg{1,x} = data_all_avg{1,x}(:, columnsToKeep);
end


%% Load the peak-definition table (collapsed-localizer output)
% This Excel file is produced upstream on the grand-grand-average across
% all participants and diagnostic groups. For each peak it specifies the
% stimulation area, the polarity-encoded peak name, the time window, and
% the four-electrode ROI. These values are fixed across subjects.
peak_excel = readtable( ...
    'F:\gbertazz\BERTANAS\giacomo\disco_F\GR2016_20230119\04_analisi_stat\analisi_TEP_peaks\peak_analysis_202406\peak_extraction_grand_grand_average.xlsx', ...
    'Sheet', 'TEP_test', ...
    'PreserveVariableNames', 1);
peak_excel_str = table2struct(peak_excel);


%% Flatten dataset cell-of-cells into a single 5-row metadata cell array
% After this step:
%   row 1: FieldTrip timelock struct
%   row 2: dataset name (BIDS-style; no group info)
%   row 3: stimulation area (LF/RF/LP/RP)
%   row 4: subject ID
%   row 5: session
data_all_avg_sub = [data_all_avg{1, 1} data_all_avg{1, 2}  data_all_avg{1, 3}  data_all_avg{1, 4}];
for y = 1:length(data_all_avg_sub)
    data_all_avg_sub{2,y} = data_all_avg_sub{1,y}.setname;
    data_all_avg_sub{3,y} = char(extractBetween(data_all_avg_sub{1,y}.setname, 'rest_', '_eeg'));
    data_all_avg_sub{4,y} = char(extractBetween(data_all_avg_sub{1,y}.setname, 'sub-',  '_ses'));
    data_all_avg_sub{5,y} = char(extractBetween(data_all_avg_sub{1,y}.setname, 'ses-',  '_task'));
end


%% Run mode
% overwright = 0 -> resume from last saved state (skip already-done
%                   subject x peak combinations).
% overwright = 1 -> redo all peaks from scratch.
overwright = 0;


%% Resume logic
% If a partial output struct already exists, load it and continue from the
% next entry. Otherwise initialize an empty struct.
if exist([savepath '\first_peak_pos_grand_gavg_avg_TP9_TP10_test.mat'], 'file') == 2
    load([savepath '\first_peak_pos_grand_gavg_avg_TP9_TP10_test'])
    global_count = size(first_peak_pos_grand_gavg_avg_TP9_TP10_test, 2) + 1;
else
    global_count = 1;
    first_peak_pos_grand_gavg_avg_TP9_TP10_test(1).name = ''; % dummy entry
end


%% Main loop: iterate over datasets (subject x stimulation site)
for dataset_count = 1:size(data_all_avg_sub, 2)


    %% Select the peaks relevant to the current stimulation area
    % Each stimulation site (LF, RF, LP, RP) is associated with a fixed
    % set of peaks defined in the peak-definition table.
    current_area       = data_all_avg_sub(3, dataset_count);
    peak_excel_str_sel = peak_excel_str(contains({peak_excel_str.area}, current_area));


    %% Inner loop: iterate over peaks for the current dataset
    for peak_count = 1:length(peak_excel_str_sel)

        % --- Retrieve fixed parameters for this peak (from the collapsed
        %     localizer, identical across subjects):
        peak_lat     = [peak_excel_str_sel(peak_count).start ...
                        peak_excel_str_sel(peak_count).xEnd];   % time window
        peak_elec    = split(peak_excel_str_sel(peak_count).elec); % ROI electrodes
        current_peak = peak_excel_str_sel(peak_count).peakName;   % e.g. 'N20F_amp_L_DLPFC'

        % --- Skip if this (subject x peak) has already been processed
        %     (resume mode).
        if  overwright == 0 ...
                && any(strcmp(data_all_avg_sub(2, dataset_count), {first_peak_pos_grand_gavg_avg_TP9_TP10_test.name})) ...
                && any(strcmp(current_peak, {first_peak_pos_grand_gavg_avg_TP9_TP10_test(strcmp(data_all_avg_sub(2, dataset_count), {first_peak_pos_grand_gavg_avg_TP9_TP10_test.name})).peakName}))

            warning(['******  dataset ' current_peak '-' char(current_area) '-' data_all_avg_sub{2,dataset_count} ' already done, skipping to the next ******']);

        else


            %% Select the data slice for this peak
            % Restrict the trace to a generous display window
            % (-10 ms to 300 ms) on the predefined ROI channels.
            cfg          = [];
            cfg.latency  = [-0.010 0.300];
            cfg.channel  = peak_elec;
            dummy_sel_data = ft_selectdata(cfg, data_all_avg_sub{1, dataset_count});

            % Time vector expressed in decimal milliseconds for exact
            % integer indexing of the search window.
            time_epoch = (round(dummy_sel_data.time * 10000));


            %% Locate the search-window indices within the time vector
            TEP                  = dummy_sel_data.avg;  % raw (non-rectified) TEP
            interval_start_index = find(time_epoch == round(peak_lat(1) * 10000), 1);
            interval_end_index   = find(time_epoch == round(peak_lat(2) * 10000), 1);


            %% Automatic candidate detection with findpeaks
            % For each ROI electrode, scan the predefined window for
            % local maxima above a minimum prominence (0.05). The largest
            % candidate across all ROI electrodes is selected as the
            % automatic pick. For negative peaks the signal is inverted.
            PKS_first_ch  = nan(50, 100);  % preallocated buffers
            LOCS_first_ch = nan(50, 100);

            if startsWith(current_peak, 'P')
                % -------------------- Positive peaks --------------------
                for find_peak_count = 1:length(dummy_sel_data.label)
                    [find_peak_amp, find_peak_lat] = ...
                        findpeaks(dummy_sel_data.avg(find_peak_count, interval_start_index:interval_end_index), ... % Y
                                  peak_lat(1):0.0002:peak_lat(2), ...                                                % X
                                  'MinPeakProminence', 0.05, ...
                                  'Annotate', 'peaks');
                    PKS_first_ch(find_peak_count,  1:length(find_peak_amp)) = find_peak_amp;
                    LOCS_first_ch(find_peak_count, 1:length(find_peak_lat)) = find_peak_lat;
                end

                % Choose the largest candidate across electrodes
                PKS_first_ch_max    = max(max(PKS_first_ch));
                [chan_index, peak_index] = find(PKS_first_ch == PKS_first_ch_max);
                LOCS_first_ch_max   = LOCS_first_ch(chan_index, peak_index);

                % Placeholder if no candidate is found
                if isempty(PKS_first_ch_max) || isnan(PKS_first_ch_max)
                    PKS_first_ch_max  = 0.000000000000001;
                    LOCS_first_ch_max = -0.010;
                end

            elseif startsWith(current_peak, 'N')
                % -------------------- Negative peaks --------------------
                for find_peak_count = 1:length(dummy_sel_data.label)
                    [find_peak_amp, find_peak_lat] = ...
                        findpeaks(-(dummy_sel_data.avg(find_peak_count, interval_start_index:interval_end_index)), ... % Y inverted
                                  peak_lat(1):0.0002:peak_lat(2), ...                                                   % X
                                  'MinPeakProminence', 0.05, ...
                                  'Annotate', 'peaks');
                    PKS_first_ch(find_peak_count,  1:length(find_peak_amp)) = find_peak_amp;
                    LOCS_first_ch(find_peak_count, 1:length(find_peak_lat)) = find_peak_lat;
                end

                % Choose the largest candidate, then restore sign
                PKS_first_ch_max    = max(max(PKS_first_ch));
                [chan_index, peak_index] = find(PKS_first_ch == PKS_first_ch_max);
                LOCS_first_ch_max   = LOCS_first_ch(chan_index, peak_index);
                PKS_first_ch_max    = -PKS_first_ch_max;  % undo earlier inversion

                % Placeholder if no candidate is found
                if isempty(PKS_first_ch_max) || isnan(PKS_first_ch_max)
                    PKS_first_ch_max  = 0.000000000000001;
                    LOCS_first_ch_max = -0.010;
                end
            end


            %% Plot the trace, highlight the search window and the pick
            % The figure is maximized for visual inspection. The window
            % boundaries are drawn as vertical lines, and the automatic
            % pick is marked with a green star.
            h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

            plot(dummy_sel_data.time, dummy_sel_data.avg)  % plot TEP
            hold on

            set(0, 'DefaulttextInterpreter', 'none')
            xline(peak_lat(1))  % search-window start
            xline(peak_lat(2))  % search-window end

            % Marker at the automatic pick
            plot(LOCS_first_ch_max, PKS_first_ch_max, 'g*', 'LineWidth', 2)

            legend(dummy_sel_data.label);
            title([peak_excel_str_sel(peak_count).peakName ' - ' ...
                   {strtok(data_all_avg_sub{1, dataset_count}.setname, '.')}]);
            subtitle('press 1 when the peak is selected. if there is not peak, select the highest point in the time window and press 0');

            hold off

            % Enable data-cursor mode so the operator can click on the
            % trace to override the automatic pick.
            datacursormode on;
            dcm_obj = datacursormode(h);


            %% Human accept / reject / override
            % Wait for the operator's key press:
            %   '1' -> accept (with or without cursor override)
            %   '0' -> reject (no peak); latency is emptied
            % Any other key is ignored and the loop continues.
            key_check = 0;
            while key_check == 0
                c_info = getCursorInfo(dcm_obj);  % cursor info, if any
                warning('press 1 or 0')
                waitforbuttonpress
                if str2double(get(gcf, 'CurrentCharacter')) == 1
                    break                          % accept
                elseif str2double(get(gcf, 'CurrentCharacter')) == 0
                    LOCS_first_ch_max = [];        % reject -> empty latency
                    break
                else
                    % ignore other keys
                end
            end


            %% Store the result, depending on operator action
            % Three cases:
            %   (a) no cursor was used -> accept the automatic pick;
            %   (b) cursor was used but empty -> still accept automatic pick;
            %   (c) cursor was used and non-empty -> use the manually
            %       clicked coordinates as the override.
            if exist('c_info', 'var') == 0
                % Case (a): operator accepted the automatic pick
                first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).name     = strtok(data_all_avg_sub{1, dataset_count}.setname, '.');
                first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).ampli    = PKS_first_ch_max;
                first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).latency  = LOCS_first_ch_max;
                first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).channel  = dummy_sel_data.label{chan_index};
                first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).peakName = peak_excel_str_sel(peak_count).peakName

            elseif exist('c_info', 'var') == 1 && isempty(c_info) == 1
                % Case (b): cursor was queried but not placed; accept automatic pick
                first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).name     = strtok(data_all_avg_sub{1, dataset_count}.setname, '.');
                first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).ampli    = PKS_first_ch_max;
                first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).latency  = LOCS_first_ch_max;
                first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).channel  = dummy_sel_data.label{chan_index};
                first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).peakName = peak_excel_str_sel(peak_count).peakName

            elseif exist('c_info', 'var') == 1 && isempty(c_info) == 0
                % Case (c): operator manually overrode the automatic pick
                first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).name     = strtok(data_all_avg_sub{1, dataset_count}.setname, '.');
                first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).ampli    = c_info.Position(2);
                first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).latency  = c_info.Position(1);
                first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).channel  = c_info.Target.DisplayName;
                first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).peakName = peak_excel_str_sel(peak_count).peakName
            end

            clear c_info  % reset cursor info for next iteration


            %% Save the inspection figure (auditability)
            savefig([savepath '/first_peak_pos_grand_gavg_avg_TP9_TP10_test_' ...
                     strtok(data_all_avg_sub{1, dataset_count}.setname, '.') '_' ...
                     peak_excel_str_sel(peak_count).peakName]);
            saveas(h, ...
                  [savepath '/first_peak_pos_grand_gavg_avg_TP9_TP10_test_' ...
                   strtok(data_all_avg_sub{1, dataset_count}.setname, '.') '_' ...
                   peak_excel_str_sel(peak_count).peakName], ...
                  'tif');
            close all


            %% Save the running results after every accept/reject
            % This makes the procedure resumable in case of interruption.
            save([savepath '/first_peak_pos_grand_gavg_avg_TP9_TP10_test'], ...
                 'first_peak_pos_grand_gavg_avg_TP9_TP10_test');

            close all
            global_count = global_count + 1;
            clc
        end

    end
end


%% =====================================================================
%  RESHAPE: long-format struct -> wide-format table (subjects x measures)
%  =====================================================================
% Convert the per-event struct (one entry per subject x peak) into a wide
% table with one row per subject and one set of columns
% (amplitude, latency, electrode) per (peak x stimulation area).

% Reload the final struct from disk (defensive: same file just written
% above, but loading it makes this block runnable independently).
load('F:\gbertazz\BERTANAS\giacomo\disco_F\GR2016_20230119\04_analisi_stat\analisi_TEP_peaks\peak_analysis_202501\peak_extraction\first_peak_pos_grand_gavg_avg_TP9_TP10_test.mat');

% Extract unique subject IDs from the dataset names
num_entries  = length(first_peak_pos_grand_gavg_avg_TP9_TP10_test);
subject_ids  = cell(num_entries, 1);

for i = 1:num_entries
    name_parts       = strsplit(first_peak_pos_grand_gavg_avg_TP9_TP10_test(i).name, '_');
    sub_num          = strsplit(name_parts{1}, '-');
    sub_num          = sub_num{2};
    subject_ids{i}   = [sub_num '_' name_parts{2}];
end
unique_subjects = sort(unique(subject_ids));  % sorted for consistent column order

% Extract unique peak names
all_peak_names = {first_peak_pos_grand_gavg_avg_TP9_TP10_test.peakName};
unique_peaks   = unique(all_peak_names);

% Split each peak name into base peak (e.g. 'N20F') and stimulation area
% (e.g. 'L_DLPFC') to build column names.
base_peaks = cell(size(unique_peaks));
stim_areas = cell(size(unique_peaks));
for i = 1:length(unique_peaks)
    peak_parts    = strsplit(unique_peaks{i}, '_');
    base_peaks{i} = peak_parts{1};
    stim_areas{i} = [peak_parts{end-1} '_' peak_parts{end}];
end
unique_base_peaks = sort(unique(base_peaks));
unique_stim_areas = sort(unique(stim_areas));

% Initialize the wide-format struct
new_struct            = struct();
new_struct.subject_id = unique_subjects;

% Map subject IDs -> row indices (faster than repeated strcmp)
subject_to_idx = containers.Map(unique_subjects, 1:length(unique_subjects));

% Create three columns per (peak x area): amplitude, latency, electrode
for i = 1:length(unique_base_peaks)
    for j = 1:length(unique_stim_areas)
        base_peak = unique_base_peaks{i};
        stim_area = unique_stim_areas{j};

        field_amp  = ['T0_' base_peak '_amp_'  stim_area];
        field_lat  = ['T0_' base_peak '_lat_'  stim_area];
        field_elec = ['T0_' base_peak '_elec_' stim_area];

        new_struct.(field_amp)  = nan(length(unique_subjects), 1);
        new_struct.(field_lat)  = nan(length(unique_subjects), 1);
        new_struct.(field_elec) = cell(length(unique_subjects), 1);
        [new_struct.(field_elec){:}] = deal('');  % default empty string
    end
end

% Populate the wide-format struct from the long-format records
for idx = 1:num_entries
    peak_info = first_peak_pos_grand_gavg_avg_TP9_TP10_test(idx);

    % Rebuild the subject ID
    name_parts  = strsplit(peak_info.name, '_');
    sub_num     = strsplit(name_parts{1}, '-');
    subject_id  = [sub_num{2} '_' name_parts{2}];
    subject_idx = subject_to_idx(subject_id);

    % Decompose the peak name into base peak and stimulation area
    peak_parts = strsplit(peak_info.peakName, '_');
    base_peak  = peak_parts{1};
    stim_area  = [peak_parts{end-1} '_' peak_parts{end}];

    field_amp  = ['T0_' base_peak '_amp_'  stim_area];
    field_lat  = ['T0_' base_peak '_lat_'  stim_area];
    field_elec = ['T0_' base_peak '_elec_' stim_area];

    % Write values (only if present in the struct)
    if ~isempty(peak_info.ampli)
        new_struct.(field_amp)(subject_idx) = peak_info.ampli;
    end
    if ~isempty(peak_info.latency)
        new_struct.(field_lat)(subject_idx) = peak_info.latency;
    end
    if ~isempty(peak_info.channel)
        new_struct.(field_elec){subject_idx} = peak_info.channel;
    end
end

% Convert to a MATLAB table for export
T        = struct2table(new_struct);
colNames = T.Properties.VariableNames;

% -------------------------------------------------------------------------
% Optional polarity-consistency filter (DISABLED in the final pipeline).
% Originally available to coerce wrong-polarity picks to NaN
% (N20 components forced to be negative, P20 components forced to be
% positive). Disabled because the manual accept/reject step already
% catches polarity-inconsistent picks. Kept here as a transparent record
% of all options considered.
% -------------------------------------------------------------------------
% % Process N20 peaks
% for i = 1:length(colNames)
%     currentCol = colNames{i};
%
%     % Check if it's an N20 amplitude column
%     if contains(currentCol, 'N20') && contains(currentCol, '_amp_')
%         % Find corresponding latency and channel columns
%         latCol  = strrep(currentCol, '_amp_', '_lat_');
%         elecCol = strrep(currentCol, '_amp_', '_elec_');
%
%         % Find indices where amplitude is > 0
%         idx = T.(currentCol) > 0;
%
%         % Set values to NaN/empty where condition is met
%         T.(currentCol)(idx) = NaN;
%         T.(latCol)(idx)     = NaN;
%         T.(elecCol)(idx)    = {''};
%     end
%
%     % Check if it's a P20 amplitude column
%     if contains(currentCol, 'P20') && contains(currentCol, '_amp_')
%         % Find corresponding latency and channel columns
%         latCol  = strrep(currentCol, '_amp_', '_lat_');
%         elecCol = strrep(currentCol, '_amp_', '_elec_');
%
%         % Find indices where amplitude is < 0
%         idx = T.(currentCol) < 0;
%
%         % Set values to NaN/empty where condition is met
%         T.(currentCol)(idx) = NaN;
%         T.(latCol)(idx)     = NaN;
%         T.(elecCol)(idx)    = {''};
%     end
% end


%% Save final outputs
% .mat -> long-format struct (one entry per subject x peak)
save([savepath '/first_peak_pos_grand_gavg_avg_TP9_TP10_test'], ...
     'first_peak_pos_grand_gavg_avg_TP9_TP10_test');

% .csv -> wide-format table for downstream statistical analysis in R
writetable(T, [savepath '/first_peak_pos_grand_gavg_avg_TP9_TP10_test.csv']);
