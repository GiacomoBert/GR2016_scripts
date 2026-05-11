% =========================================================================
% peak_extraction.m
% -------------------------------------------------------------------------
% Semi-automatic extraction of early TMS-evoked potential (TEP) peaks,
% following a collapsed-localizer strategy seeded on the grand-grand-average.
%
% Companion code for:
%   Bertazzoli G., Canu E., Bagattini C., et al.
%   "Network-targeted neurophysiological biomarkers of dysconnectivity and
%   cognitive decline in Alzheimer's disease."
%   https://doi.org/10.5281/zenodo.17162745
%
% Author : Giacomo Bertazzoli (gbertazz@bidmc.harvard.edu)
% License: CC-BY-4.0
%
% Procedure:
%   1. Time window + ROI electrodes for each peak are defined on the
%      grand-grand-average (all subjects, all groups) and stored in
%      `peak_extraction_grand_grand_average.xlsx`. They are fixed across
%      subjects and independent of diagnostic group.
%   2. Per subject and peak: `findpeaks` runs on each ROI electrode within
%      the predefined window (MinPeakProminence = 0.05; signal inverted
%      for negative peaks). Largest-prominence candidate is the automatic
%      pick.
%   3. Operator decides: '1' accept, '0' reject (no peak), or click on
%      the trace then '1' to override. Titles/filenames carry only
%      anonymized subject IDs (no group info).
%   4. Saved after every decision (resumable). Final long-format struct
%      is reshaped into a wide CSV.
%
% NOTE: the canonical search windows and ROI electrodes used in the
% published analysis are those reported in the paper (Table 2 / Methods)
% and on Zenodo (DOI above). The local spreadsheet should be edited to
% match these before running.
%
% Dependencies: MATLAB R2020a+, EEGLAB, FieldTrip, natsortfiles.
% =========================================================================


%% Setup
clear all
close all
clc
main_dir = extractBefore(mfilename('fullpath'), ['\' mfilename]);


%% Plug-ins (adjust paths to your local installation)
addpath('E:\BrightFocus_PreclinicalAD\01_preprocessing\plugins\natsortfiles');
addpath('E:\BrightFocus_PreclinicalAD\01_preprocessing\plugins\eeglab2024.2');
addpath('E:\BrightFocus_PreclinicalAD\01_preprocessing\plugins\fieldtrip-20240731');


%% Output folder
savepath = [main_dir '/peak_extraction'];
mkdir(savepath);


%% Init EEGLAB / FieldTrip
eeglab nogui
ft_defaults


%% Paths to source data (for reference; not read by this script)
raw_folder  = 'F:\gbertazz\BERTANAS\archive\raw_data';
data_parent = 'F:\gbertazz\BERTANAS\archive\cleaned_data\derivatives';


%% Load subject-level grand averages
% Cell array, one cell per stimulation site (LF/RF/LP/RP); each cell
% contains FieldTrip timelock structs (one per subject x session).
load('F:\gbertazz\BERTANAS\giacomo\disco_F\GR2016_20230119\04_analisi_stat\grand_averages\gavg\grand_gavg_avg_2024\data_all_avg.mat');


%% Restrict to subjects in the master list
filename   = 'H:\GR2016_peer\TEP_DTI_measures_total_20240315.xlsx';
sheet      = 'Sheet1';
data       = readtable(filename, 'sheet', sheet);
subj       = table2cell(data(:, 'ID_subj_NT'));
subj       = extractBefore(subj, 7);                 % strip group letter; keep anonymous ID
filterStrings = subj;

for x = 1:length(data_all_avg)
    for y = 1:length(data_all_avg{1,x})
        if ~contains(data_all_avg{1,x}{y}.setname, filterStrings)
            data_all_avg{1,x}{y} = [];
        end
    end
    columnsToKeep = sum(cellfun(@isempty, data_all_avg{1,x}), 1) ~= size(data_all_avg{1,x}, 1);
    data_all_avg{1,x} = data_all_avg{1,x}(:, columnsToKeep);
end


%% Load peak-definition table (collapsed-localizer output)
% Columns: area, peakName, start, xEnd, elec. Fixed across subjects.
peak_excel = readtable( ...
    'F:\gbertazz\BERTANAS\giacomo\disco_F\GR2016_20230119\04_analisi_stat\analisi_TEP_peaks\peak_analysis_202406\peak_extraction_grand_grand_average.xlsx', ...
    'Sheet', 'TEP_test', 'PreserveVariableNames', 1);
peak_excel_str = table2struct(peak_excel);


%% Flatten datasets and extract metadata
% After this: row 1 = struct, row 2 = setname, row 3 = stim area,
% row 4 = subject ID, row 5 = session.
data_all_avg_sub = [data_all_avg{1,1} data_all_avg{1,2} data_all_avg{1,3} data_all_avg{1,4}];
for y = 1:length(data_all_avg_sub)
    data_all_avg_sub{2,y} = data_all_avg_sub{1,y}.setname;
    data_all_avg_sub{3,y} = char(extractBetween(data_all_avg_sub{1,y}.setname, 'rest_', '_eeg'));
    data_all_avg_sub{4,y} = char(extractBetween(data_all_avg_sub{1,y}.setname, 'sub-',  '_ses'));
    data_all_avg_sub{5,y} = char(extractBetween(data_all_avg_sub{1,y}.setname, 'ses-',  '_task'));
end


%% Run mode (0 = resume, 1 = redo all)
overwright = 0;


%% Resume from saved state if present
if exist([savepath '\first_peak_pos_grand_gavg_avg_TP9_TP10_test.mat'], 'file') == 2
    load([savepath '\first_peak_pos_grand_gavg_avg_TP9_TP10_test'])
    global_count = size(first_peak_pos_grand_gavg_avg_TP9_TP10_test, 2) + 1;
else
    global_count = 1;
    first_peak_pos_grand_gavg_avg_TP9_TP10_test(1).name = '';
end


%% Main loop: subject x stimulation site
for dataset_count = 1:size(data_all_avg_sub, 2)

    current_area       = data_all_avg_sub(3, dataset_count);
    peak_excel_str_sel = peak_excel_str(contains({peak_excel_str.area}, current_area));

    for peak_count = 1:length(peak_excel_str_sel)

        % Fixed parameters for this peak
        peak_lat     = [peak_excel_str_sel(peak_count).start peak_excel_str_sel(peak_count).xEnd];
        peak_elec    = split(peak_excel_str_sel(peak_count).elec);
        current_peak = peak_excel_str_sel(peak_count).peakName;

        % Skip if already done (resume mode)
        if  overwright == 0 ...
                && any(strcmp(data_all_avg_sub(2, dataset_count), {first_peak_pos_grand_gavg_avg_TP9_TP10_test.name})) ...
                && any(strcmp(current_peak, {first_peak_pos_grand_gavg_avg_TP9_TP10_test(strcmp(data_all_avg_sub(2, dataset_count), {first_peak_pos_grand_gavg_avg_TP9_TP10_test.name})).peakName}))
            warning(['******  ' current_peak '-' char(current_area) '-' data_all_avg_sub{2,dataset_count} ' already done ******']);
            continue
        end


        %% Slice data: -10 to 300 ms on ROI channels
        cfg = []; cfg.latency = [-0.010 0.300]; cfg.channel = peak_elec;
        dummy_sel_data = ft_selectdata(cfg, data_all_avg_sub{1, dataset_count});

        time_epoch = round(dummy_sel_data.time * 10000);
        interval_start_index = find(time_epoch == round(peak_lat(1) * 10000), 1);
        interval_end_index   = find(time_epoch == round(peak_lat(2) * 10000), 1);


        %% Automatic candidate detection (findpeaks per ROI electrode)
        PKS_first_ch  = nan(50, 100);
        LOCS_first_ch = nan(50, 100);

        if startsWith(current_peak, 'P')
            % Positive peaks
            for find_peak_count = 1:length(dummy_sel_data.label)
                [find_peak_amp, find_peak_lat] = findpeaks( ...
                    dummy_sel_data.avg(find_peak_count, interval_start_index:interval_end_index), ...
                    peak_lat(1):0.0002:peak_lat(2), 'MinPeakProminence', 0.05, 'Annotate', 'peaks');
                PKS_first_ch(find_peak_count,  1:length(find_peak_amp)) = find_peak_amp;
                LOCS_first_ch(find_peak_count, 1:length(find_peak_lat)) = find_peak_lat;
            end
            PKS_first_ch_max = max(max(PKS_first_ch));
            [chan_index, peak_index] = find(PKS_first_ch == PKS_first_ch_max);
            LOCS_first_ch_max = LOCS_first_ch(chan_index, peak_index);
            if isempty(PKS_first_ch_max) || isnan(PKS_first_ch_max)
                PKS_first_ch_max = 0.000000000000001; LOCS_first_ch_max = -0.010;
            end

        elseif startsWith(current_peak, 'N')
            % Negative peaks: invert signal, run findpeaks, restore sign
            for find_peak_count = 1:length(dummy_sel_data.label)
                [find_peak_amp, find_peak_lat] = findpeaks( ...
                    -(dummy_sel_data.avg(find_peak_count, interval_start_index:interval_end_index)), ...
                    peak_lat(1):0.0002:peak_lat(2), 'MinPeakProminence', 0.05, 'Annotate', 'peaks');
                PKS_first_ch(find_peak_count,  1:length(find_peak_amp)) = find_peak_amp;
                LOCS_first_ch(find_peak_count, 1:length(find_peak_lat)) = find_peak_lat;
            end
            PKS_first_ch_max = max(max(PKS_first_ch));
            [chan_index, peak_index] = find(PKS_first_ch == PKS_first_ch_max);
            LOCS_first_ch_max = LOCS_first_ch(chan_index, peak_index);
            PKS_first_ch_max  = -PKS_first_ch_max;
            if isempty(PKS_first_ch_max) || isnan(PKS_first_ch_max)
                PKS_first_ch_max = 0.000000000000001; LOCS_first_ch_max = -0.010;
            end
        end


        %% Plot trace + window + automatic pick
        h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        plot(dummy_sel_data.time, dummy_sel_data.avg); hold on
        set(0, 'DefaulttextInterpreter', 'none')
        xline(peak_lat(1)); xline(peak_lat(2));
        plot(LOCS_first_ch_max, PKS_first_ch_max, 'g*', 'LineWidth', 2)
        legend(dummy_sel_data.label);
        title([peak_excel_str_sel(peak_count).peakName ' - ' {strtok(data_all_avg_sub{1, dataset_count}.setname, '.')}]);
        subtitle('press 1 to accept; click + 1 to override; 0 to reject');
        hold off

        datacursormode on;
        dcm_obj = datacursormode(h);


        %% Operator: 1 = accept, 0 = reject, click + 1 = override
        key_check = 0;
        while key_check == 0
            c_info = getCursorInfo(dcm_obj);
            warning('press 1 or 0')
            waitforbuttonpress
            if str2double(get(gcf, 'CurrentCharacter')) == 1
                break
            elseif str2double(get(gcf, 'CurrentCharacter')) == 0
                LOCS_first_ch_max = [];
                break
            end
        end


        %% Store result
        if exist('c_info', 'var') == 0 || isempty(c_info)
            % Accept automatic pick
            first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).name     = strtok(data_all_avg_sub{1, dataset_count}.setname, '.');
            first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).ampli    = PKS_first_ch_max;
            first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).latency  = LOCS_first_ch_max;
            first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).channel  = dummy_sel_data.label{chan_index};
            first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).peakName = peak_excel_str_sel(peak_count).peakName;
        else
            % Manual override
            first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).name     = strtok(data_all_avg_sub{1, dataset_count}.setname, '.');
            first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).ampli    = c_info.Position(2);
            first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).latency  = c_info.Position(1);
            first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).channel  = c_info.Target.DisplayName;
            first_peak_pos_grand_gavg_avg_TP9_TP10_test(global_count).peakName = peak_excel_str_sel(peak_count).peakName;
        end
        clear c_info


        %% Save figure and running results
        savefig([savepath '/first_peak_pos_grand_gavg_avg_TP9_TP10_test_' ...
                 strtok(data_all_avg_sub{1, dataset_count}.setname, '.') '_' ...
                 peak_excel_str_sel(peak_count).peakName]);
        saveas(h, [savepath '/first_peak_pos_grand_gavg_avg_TP9_TP10_test_' ...
                   strtok(data_all_avg_sub{1, dataset_count}.setname, '.') '_' ...
                   peak_excel_str_sel(peak_count).peakName], 'tif');
        close all
        save([savepath '/first_peak_pos_grand_gavg_avg_TP9_TP10_test'], ...
             'first_peak_pos_grand_gavg_avg_TP9_TP10_test');
        global_count = global_count + 1;
        clc

    end
end


%% =====================================================================
%  Reshape long-format struct -> wide-format CSV (subject x measures)
%  =====================================================================
load('F:\gbertazz\BERTANAS\giacomo\disco_F\GR2016_20230119\04_analisi_stat\analisi_TEP_peaks\peak_analysis_202501\peak_extraction\first_peak_pos_grand_gavg_avg_TP9_TP10_test.mat');

% Unique subject IDs
num_entries = length(first_peak_pos_grand_gavg_avg_TP9_TP10_test);
subject_ids = cell(num_entries, 1);
for i = 1:num_entries
    name_parts     = strsplit(first_peak_pos_grand_gavg_avg_TP9_TP10_test(i).name, '_');
    sub_num        = strsplit(name_parts{1}, '-');
    subject_ids{i} = [sub_num{2} '_' name_parts{2}];
end
unique_subjects = sort(unique(subject_ids));

% Unique peaks (split into base peak + stim area)
all_peak_names = {first_peak_pos_grand_gavg_avg_TP9_TP10_test.peakName};
unique_peaks   = unique(all_peak_names);
base_peaks = cell(size(unique_peaks));
stim_areas = cell(size(unique_peaks));
for i = 1:length(unique_peaks)
    peak_parts    = strsplit(unique_peaks{i}, '_');
    base_peaks{i} = peak_parts{1};
    stim_areas{i} = [peak_parts{end-1} '_' peak_parts{end}];
end
unique_base_peaks = sort(unique(base_peaks));
unique_stim_areas = sort(unique(stim_areas));

% Init wide struct (amp / lat / elec per peak x area)
new_struct            = struct();
new_struct.subject_id = unique_subjects;
subject_to_idx        = containers.Map(unique_subjects, 1:length(unique_subjects));
for i = 1:length(unique_base_peaks)
    for j = 1:length(unique_stim_areas)
        bp = unique_base_peaks{i}; sa = unique_stim_areas{j};
        new_struct.(['T0_' bp '_amp_'  sa]) = nan(length(unique_subjects), 1);
        new_struct.(['T0_' bp '_lat_'  sa]) = nan(length(unique_subjects), 1);
        new_struct.(['T0_' bp '_elec_' sa]) = cell(length(unique_subjects), 1);
        [new_struct.(['T0_' bp '_elec_' sa]){:}] = deal('');
    end
end

% Fill in values
for idx = 1:num_entries
    peak_info   = first_peak_pos_grand_gavg_avg_TP9_TP10_test(idx);
    name_parts  = strsplit(peak_info.name, '_');
    sub_num     = strsplit(name_parts{1}, '-');
    subject_idx = subject_to_idx([sub_num{2} '_' name_parts{2}]);

    peak_parts = strsplit(peak_info.peakName, '_');
    bp = peak_parts{1};
    sa = [peak_parts{end-1} '_' peak_parts{end}];

    if ~isempty(peak_info.ampli),   new_struct.(['T0_' bp '_amp_'  sa])(subject_idx) = peak_info.ampli;   end
    if ~isempty(peak_info.latency), new_struct.(['T0_' bp '_lat_'  sa])(subject_idx) = peak_info.latency; end
    if ~isempty(peak_info.channel), new_struct.(['T0_' bp '_elec_' sa]){subject_idx} = peak_info.channel; end
end


%% Save outputs
T = struct2table(new_struct);
save([savepath '/first_peak_pos_grand_gavg_avg_TP9_TP10_test'], 'first_peak_pos_grand_gavg_avg_TP9_TP10_test');
writetable(T, [savepath '/first_peak_pos_grand_gavg_avg_TP9_TP10_test.csv']);
