%%
% dsi_events_analysis.m
%
% Eddy Albarran (2017)
%
% This script analyzes experiments probing the magnitude of depolarization
% induced suppression of inhibition (DSi). The data should be of the form
% [pre, depol, post] - a gap-free recording of a neuron voltage clamped at
% a baseline potential, depolarized for some time, and recorded afterward.
%
% This code was used for analysis of Dong et al. (2021)
%

clear all; close all; clear_all_figs;
warning('off', 'MATLAB:unknownObjectNowStruct');
warning('off', 'MATLAB:polyfit:RepeatedPoints');

% Optional flags for analysis

PLOT_INDIVIDUAL_TRACES = 0;
PLOT_INDIVIDUAL_CELLS = 0;
PLOT_SUMMARY = 1;
DEBUG = 0;

% Signal parameters

sample_rate = 10; % 10kHz (every 10 data points is 1ms)

test1_range = [1, 300*sample_rate];
test2_range = [119400*sample_rate+1 , 1200000 - (300*sample_rate)];
pre_range = [50000*sample_rate+1 , 60000*sample_rate];
pre_length = pre_range(2)-pre_range(1)+1;
depol_range = [60000*sample_rate, 65000*sample_rate];
post_range = [65050*sample_rate+1, 95050*sample_rate];
post_length = post_range(2)-post_range(1)+1;

% Signal processing paramters

poly_order = 40; % order of polynomial fit to data for baseline subtraction
std_thresh = 1.8; % threshold for event detection (dse = 2.2, dsi = 1.8)
bin_size = 20000; % 5000 -> 500ms bins
norm_baseline_win = 100000; % num of data points to avg for norm baseline
baseline_win_binned = norm_baseline_win/bin_size;

% Plotting parameters

%data_color = [0 0 0]; 
data_color = [.8 .2 .2];
X_lims = [pre_range(2)-norm_baseline_win+1 post_range(2)];
Y_lim_freq = 12; % hz (dse = 9, dsi = 12)
Y_lim_amp = 90; % pA (dse = 50, dsi = 90)
Y_lim_charge = 12; % pC (dse = 2, dsi = 12)
Y_lim_norm = 200;
pre_range_binned = pre_range(1) : bin_size : pre_range(2);
post_range_binned = post_range(1) : bin_size : post_range(2);
baseline_plot_pts_binned = (pre_range(2)-norm_baseline_win + bin_size/2 + 1) : bin_size : pre_range(2);
post_plot_pts_binned = post_range_binned + bin_size/2;
bar_bin_size = 2;
post_end_t = 280000; % timepoint to compare return to baseline
tick_range_min = 400001;
tick_range_max = 1100001;
tick_range_step = 100000;

% Where is the data?

proj_path = [''];
data_path = [proj_path, '/data'];

data_list_path = strcat(proj_path, '');
[~, ~, data_table] = xlsread(data_list_path);

% Variables for indexing info from data table
ANIMAL = 1;
DOB = 2;
DOE = 3;
FILENAME = 4;
EXCLUDE = 5;
SOFTWARE = 6;
SEX = 7;
TRACES = 8;

%%

% Exclude data
[data_table, numRecs] = exclude_data(data_table, EXCLUDE);

% Determine the number of animals
numAnimals = length(unique([data_table{:,ANIMAL}]));

% Data structs for event freq, amp, and charge (for pre and post)

cells_pre_binned_freq = zeros(numRecs, pre_length / bin_size);
cells_pre_binned_amp = zeros(numRecs, size(cells_pre_binned_freq, 2));
cells_pre_binned_charge = zeros(numRecs, size(cells_pre_binned_freq, 2));

cells_post_binned_freq = zeros(numRecs, post_length / bin_size);
cells_post_binned_amp = zeros(numRecs, size(cells_post_binned_freq, 2));
cells_post_binned_charge = zeros(numRecs, size(cells_post_binned_freq, 2));

% Iterate through each cell, averaging the values extracted for each of the
% cell's traces
for cell_i = 1:numRecs
    
    % Parse the directory and file
    
    exp_date = num2str(data_table{cell_i, DOE});
    if length(exp_date) < 6
        exp_date = strcat('0', exp_date);
    end
    
    exp_filename = data_table{cell_i, FILENAME};
    
    abf_folder = [data_path, '/', exp_date, '/', exp_filename];
    cd(abf_folder);
    
    % Allocate structures for traces of this cell
    
    trace_indeces = num2str(data_table{cell_i, TRACES});
    num_traces = length(trace_indeces);
    
    pre_binned_freq = zeros(num_traces, pre_length / bin_size);
    pre_binned_amp = zeros(num_traces, pre_length / bin_size);
    pre_binned_charge = zeros(num_traces, pre_length / bin_size);
    
    post_binned_freq = zeros(num_traces, post_length / bin_size);
    post_binned_amp = zeros(num_traces, post_length / bin_size);
    post_binned_charge = zeros(num_traces, post_length / bin_size);
    
    % For each trace: extract freq, amp, and charge
    for trace_i = 1:num_traces
        
        % Load the trace
        
        filename = ['AD0_', num2str(trace_indeces(trace_i)), '.mat'];
        load(filename);
        eval(['d = AD0_', num2str(trace_indeces(trace_i)), '.data;']);
        
        % Determine the Ra and Ih at the beginning and end of the trace
        
        test1 = d(test1_range(1):test1_range(2));
        test2 = d(test2_range(1):test2_range(2));
        
        [Ra1, ~] = test_pulse(test1);
        [Ra2, ~] = test_pulse(test2);
        
        Ih1 = mean(test1(1:500));
        Ih2 = mean(test2(1:500));
        
        % Get baseline trace (detrend drift with polynomial fit)
        pre_depol = d(pre_range(1):pre_range(2));
        dt_pre_depol = detrend_poly_fit(pre_depol, poly_order/10);
        
        % Get post-depolarization trace (remove drift with polynomial fit)
        post_depol = d(post_range(1):post_range(2));
        dt_post_depol = detrend_poly_fit(post_depol, poly_order);
        
        % Get mean and std for baseline (pre trace)
        baseline_mean = mean(pre_depol);
        baseline_std = std(pre_depol);
        
        % Filter signals (and detrend to make event detection easier)
        [B1, A1] = butter(2, 1/100, 'low');
        
        ft_pre_depol = filter(B1, A1, pre_depol);
        ft_pre_depol = detrend_poly_fit(ft_pre_depol, poly_order);
        
        ft_post_depol = filter(B1, A1, post_depol);
        ft_post_depol = detrend_poly_fit(ft_post_depol, poly_order);
        
        % Take the derivative (for event detection)
        deriv_pre_depol = diff(ft_pre_depol);
        deriv_post_depol = diff(ft_post_depol);
        
        % Get the [binned] freq, amp, and charge for the pre and post traces
        % (see detect_events.m)
        
        ft_baseline_mean = mean(ft_pre_depol);
        ft_baseline_std = std(ft_pre_depol);
        deriv_baseline_mean = mean(deriv_pre_depol);
        deriv_baseline_std = std(deriv_pre_depol);
        deriv_pre_depol = [deriv_pre_depol, deriv_baseline_mean];
        deriv_post_depol = [deriv_post_depol, deriv_baseline_mean];
        event_thresh = deriv_baseline_mean - deriv_baseline_std * std_thresh;
        
        [cur_pre_binned_freq, cur_pre_binned_amp, cur_pre_binned_charge] = detect_events(dt_pre_depol, deriv_pre_depol, bin_size, event_thresh, 'all', 'neg');
        [cur_post_binned_freq, cur_post_binned_amp, cur_post_binned_charge] = detect_events(dt_post_depol, deriv_post_depol, bin_size, event_thresh, 'all', 'neg');
        
        % Store the variables into the superstructures
        
        pre_binned_freq(trace_i, :) = cur_pre_binned_freq;
        pre_binned_amp(trace_i, :) = cur_pre_binned_amp;
        pre_binned_charge(trace_i, :) = cur_pre_binned_charge;
        
        post_binned_freq(trace_i, :) = cur_post_binned_freq;
        post_binned_amp(trace_i, :) = cur_post_binned_amp;
        post_binned_charge(trace_i, :) = cur_post_binned_charge;
        
        % If desired, plot the individual trace
        if PLOT_INDIVIDUAL_TRACES
            
            figure('units', 'normalized', 'position', [.01 .55 .98 .35]);
            
            % Plot the pre trace (detrended as well)
            subplot(2,5,1); hold on;
            plot(pre_depol);
            plot(1:pre_length, ones(1, pre_length)*baseline_mean); plot(1:pre_length, ones(1, pre_length)*(baseline_mean - baseline_std*std_thresh));
            title('Before depolarization'); ylabel('pA');
            
            subplot(2,5,6); hold on;
            plot(dt_pre_depol);
            plot(1:pre_length, ones(1, pre_length)*ft_baseline_mean); plot(1:pre_length, ones(1, pre_length)*(ft_baseline_mean - ft_baseline_std*std_thresh));
            ylim([-150 50]);
            xlabel('time (10kHz)'); ylabel('pA (detrended)');
            
            % Plot the post trace (detrended as well)
            subplot(2,5,2); hold on;
            plot(post_depol);
            plot(1:post_length, ones(1, post_length)*baseline_mean); plot(1:post_length, ones(1, post_length)*(baseline_mean - baseline_std*std_thresh));
            title('After depolarization');
            
            subplot(2,5,7); hold on;
            plot(dt_post_depol);
            plot(1:post_length, ones(1, post_length)*ft_baseline_mean); plot(1:post_length, ones(1, post_length)*(ft_baseline_mean - ft_baseline_std*std_thresh));
            ylim([-150 50]);
            xlabel('time (10kHz)');
            
            % Plot the binned frequency for this trace
            subplot(2,5,[3,8]); hold on;
            scatter(baseline_plot_pts_binned, cur_pre_binned_freq(1 + end - baseline_win_binned : end));
            scatter(post_range_binned, cur_post_binned_freq);
            xlabel('time (10kHz)'); ylabel('Frequency (hz)');
            xlim(X_lims); ylim([0 Y_lim_freq]);
            
            % Plot the binned amplitude for this trace
            subplot(2,5,[4,9]); hold on;
            scatter(baseline_plot_pts_binned, cur_pre_binned_amp(1 + end - baseline_win_binned : end));
            scatter(post_range_binned, cur_post_binned_amp);
            title(['Cell: ', num2str(cell_i), ' Trace: ', num2str(trace_i)]);
            xlabel('time (10kHz)'); ylabel('Amplitude (pA)');
            xlim(X_lims); ylim([0 Y_lim_amp]);
            
            % Plot the binned charge transfer for this trace
            subplot(2,5,[5,10]); hold on;
            scatter(baseline_plot_pts_binned, cur_pre_binned_charge(1 + end - baseline_win_binned : end));
            scatter(post_range_binned, cur_post_binned_charge);
            xlabel('time (10kHz)'); ylabel('Charge transfer (pC)');
            xlim(X_lims); ylim([0 Y_lim_charge]);
            
        end
        
    end
    
    % Average the numbers across traces and store into superstructure
    
    cells_pre_binned_freq(cell_i, :) = nanmean(pre_binned_freq, 1);
    cells_pre_binned_amp(cell_i, :) = nanmean(pre_binned_amp, 1);
    cells_pre_binned_charge(cell_i, :) = nanmean(pre_binned_charge, 1);
    
    cells_post_binned_freq(cell_i, :) = nanmean(post_binned_freq, 1);
    cells_post_binned_amp(cell_i, :) = nanmean(post_binned_amp, 1);
    cells_post_binned_charge(cell_i, :) = nanmean(post_binned_charge, 1);
    
    % If desired, plot the individual cells (averaged traces)
    if PLOT_INDIVIDUAL_CELLS
        
        figure('units', 'normalized', 'position', [.01 .1 .98 .35]);
        
        % Plot the binned frequency for this cell
        subplot(1,3,1); hold on;
        scatter(baseline_plot_pts_binned, cells_pre_binned_freq(cell_i, 1 + end - baseline_win_binned : end));
        scatter(post_range_binned, cells_post_binned_freq(cell_i, :));
        xlim(X_lims); ylim([0 Y_lim_freq]);
        xlabel('Time (10kHz)'); ylabel('Frequency (hz)');
        
        % Plot the binned amplitude for this cell
        subplot(1,3,2); hold on;
        scatter(baseline_plot_pts_binned, cells_pre_binned_amp(cell_i, 1 + end - baseline_win_binned : end));
        scatter(post_range_binned, cells_post_binned_amp(cell_i, :));
        title(['Cell: ', num2str(cell_i), ' (', num2str(num_traces), ' traces)']);
        xlim(X_lims); ylim([0 Y_lim_amp]);
        xlabel('Time (10kHz)'); ylabel('Amplitude (pA)');
        
        % Plot the binned charge transfer for this cell
        subplot(1,3,3); hold on;
        scatter(baseline_plot_pts_binned, cells_pre_binned_charge(cell_i, 1 + end - baseline_win_binned : end));
        scatter(post_range_binned, cells_post_binned_charge(cell_i, :));
        xlim(X_lims); ylim([0 Y_lim_charge]);
        xlabel('Time (10kHz)'); ylabel('Charge transfer (pC)');
        
    end
    
end

%%

% Calculate average values across cells for binned freq, amp, and charge

cells_pre_binned_freq_mean = mean(cells_pre_binned_freq, 1);
cells_pre_binned_amp_mean = mean(cells_pre_binned_amp, 1);
cells_pre_binned_charge_mean = mean(cells_pre_binned_charge, 1);

cells_post_binned_freq_mean = mean(cells_post_binned_freq, 1);
cells_post_binned_amp_mean = mean(cells_post_binned_amp, 1);
cells_post_binned_charge_mean = mean(cells_post_binned_charge, 1);

% Calculate SEM across cells for binned freq, amp, and charge

cells_pre_binned_freq_SEM = std(cells_pre_binned_freq, [], 1) ./ sqrt(sum(~isnan(cells_pre_binned_freq)));
cells_pre_binned_amp_SEM = std(cells_pre_binned_amp, [], 1) ./ sqrt(sum(~isnan(cells_pre_binned_amp)));
cells_pre_binned_charge_SEM = std(cells_pre_binned_charge, [], 1) ./ sqrt(sum(~isnan(cells_pre_binned_charge)));
    
cells_post_binned_freq_SEM = std(cells_post_binned_freq, [], 1) ./ sqrt(sum(~isnan(cells_post_binned_freq)));
cells_post_binned_amp_SEM = std(cells_post_binned_amp, [], 1) ./ sqrt(sum(~isnan(cells_post_binned_amp)));
cells_post_binned_charge_SEM = std(cells_post_binned_charge, [], 1) ./ sqrt(sum(~isnan(cells_post_binned_charge)));

%%
% Normalize the cell data (relative to average value of last 10 sec of pre)

norm_cells_pre_binned_freq = 100 * cells_pre_binned_freq ./ repmat(nanmean(cells_pre_binned_freq(:, 1+end-norm_baseline_win/bin_size:end), 2), 1, size(cells_pre_binned_freq, 2));
norm_cells_pre_binned_amp = 100 * cells_pre_binned_amp ./ repmat(nanmean(cells_pre_binned_amp(:, 1+end-norm_baseline_win/bin_size:end), 2), 1, size(cells_pre_binned_amp, 2));
norm_cells_pre_binned_charge = 100 * cells_pre_binned_charge ./ repmat(nanmean(cells_pre_binned_charge(:, 1+end-norm_baseline_win/bin_size:end), 2), 1, size(cells_pre_binned_charge, 2));
    
norm_cells_post_binned_freq = 100 * cells_post_binned_freq ./ repmat(nanmean(cells_pre_binned_freq(:, 1+end-norm_baseline_win/bin_size:end), 2), 1, size(cells_post_binned_freq, 2));
norm_cells_post_binned_amp = 100 * cells_post_binned_amp ./ repmat(nanmean(cells_pre_binned_amp(:, 1+end-norm_baseline_win/bin_size:end), 2), 1, size(cells_post_binned_amp, 2));
norm_cells_post_binned_charge = 100 * cells_post_binned_charge ./ repmat(nanmean(cells_pre_binned_charge(:, 1+end-norm_baseline_win/bin_size:end), 2), 1, size(cells_post_binned_charge, 2));

% Calculate normalized average values for binned freq, amp, and charge

norm_cells_pre_binned_freq_mean = mean(norm_cells_pre_binned_freq, 1);
norm_cells_pre_binned_amp_mean = mean(norm_cells_pre_binned_amp, 1);
norm_cells_pre_binned_charge_mean = mean(norm_cells_pre_binned_charge, 1);
    
norm_cells_post_binned_freq_mean = mean(norm_cells_post_binned_freq, 1);
norm_cells_post_binned_amp_mean = mean(norm_cells_post_binned_amp, 1);
norm_cells_post_binned_charge_mean = mean(norm_cells_post_binned_charge, 1);

% Calculate normalized SEM across cells for binned freq, amp, and charge

norm_cells_pre_binned_freq_SEM = std(norm_cells_pre_binned_freq, [], 1) ./ sqrt(sum(~isnan(norm_cells_pre_binned_freq)));
norm_cells_pre_binned_amp_SEM = std(norm_cells_pre_binned_amp, [], 1) ./ sqrt(sum(~isnan(norm_cells_pre_binned_amp)));
norm_cells_pre_binned_charge_SEM = std(norm_cells_pre_binned_charge, [], 1) ./ sqrt(sum(~isnan(norm_cells_pre_binned_charge)));
    
norm_cells_post_binned_freq_SEM = std(norm_cells_post_binned_freq, [], 1) ./ sqrt(sum(~isnan(norm_cells_post_binned_freq)));
norm_cells_post_binned_amp_SEM = std(norm_cells_post_binned_amp, [], 1) ./ sqrt(sum(~isnan(norm_cells_post_binned_amp)));
norm_cells_post_binned_charge_SEM = std(norm_cells_post_binned_charge, [], 1) ./ sqrt(sum(~isnan(norm_cells_post_binned_charge)));

% Plot the summarized data (raw values and normalized)

if PLOT_SUMMARY
    
    figure('units', 'normalized', 'position', [.01 .1 .98 .8]); hold on;
    
    % Plot the averaged binned frequency (across cells)
    
    subplot(2,3,1); hold on;
    title('Frequency'); 
    scatter(baseline_plot_pts_binned, cells_pre_binned_freq_mean(1 + end - baseline_win_binned : end), 80, data_color, 'filled');
    errbar(baseline_plot_pts_binned, cells_pre_binned_freq_mean(1 + end - baseline_win_binned : end), cells_pre_binned_freq_SEM(1 + end - baseline_win_binned : end), 'color', data_color);
    scatter(post_plot_pts_binned, cells_post_binned_freq_mean, 80, data_color, 'filled');
    errbar(post_plot_pts_binned, cells_post_binned_freq_mean, cells_post_binned_freq_SEM, 'color', data_color);
    rectangle('Position',[depol_range(1) Y_lim_freq*0.47 depol_range(2)-depol_range(1) Y_lim_freq*.03], 'facecolor', [0 0 0]);
    text(depol_range(1) + (depol_range(2)-depol_range(1))/2, Y_lim_freq*.55, 'depol', 'color', [0 0 0], 'horizontalalignment', 'center', 'fontsize', 10, 'fontweight', 'bold');
    text(1000000, 1, ['n = ', num2str(numRecs)], 'color', data_color);
    xlabel('Time (sec)'); ylabel('Frequency (hz)');
    xlim(X_lims); ylim([0 Y_lim_freq]);
    set(gca,'xtick', tick_range_min : tick_range_step : tick_range_max, 'xticklabel', {'-20','-10', '0', '10', '20', '30', '40', '50'});
    
    % Plot the averaged binned amplitude (across cells)
    
    subplot(2,3,2); hold on;
    title('Amplitude');
    scatter(baseline_plot_pts_binned, cells_pre_binned_amp_mean(1 + end - baseline_win_binned : end), 80, data_color, 'filled');
    errbar(baseline_plot_pts_binned, cells_pre_binned_amp_mean(1 + end - baseline_win_binned : end), cells_pre_binned_amp_SEM(1 + end - baseline_win_binned : end), 'color', data_color);
    scatter(post_plot_pts_binned, cells_post_binned_amp_mean, 80, data_color, 'filled');
    errbar(post_plot_pts_binned, cells_post_binned_amp_mean, cells_post_binned_amp_SEM, 'color', data_color);
    rectangle('Position',[depol_range(1) Y_lim_amp*0.47 depol_range(2)-depol_range(1) Y_lim_amp*.03], 'facecolor', [0 0 0]);
    text(depol_range(1) + (depol_range(2)-depol_range(1))/2, Y_lim_amp*.55, 'depol', 'color', [0 0 0], 'horizontalalignment', 'center', 'fontsize', 10, 'fontweight', 'bold');
    text(1000000, 8, ['n = ', num2str(numRecs)], 'color', data_color);
    xlabel('Time (sec)'); ylabel('Amplitude (pA)');
    xlim(X_lims); ylim([0 Y_lim_amp]);
    set(gca,'xtick', tick_range_min : tick_range_step : tick_range_max, 'xticklabel', {'-20','-10', '0', '10', '20', '30', '40', '50'})
    
    % Plot the averaged binned charge transfer (across cells)
    
    subplot(2,3,3); hold on;
    title('Charge Transfer'); 
    scatter(baseline_plot_pts_binned, cells_pre_binned_charge_mean(1 + end - baseline_win_binned : end), 80, data_color, 'filled');
    errbar(baseline_plot_pts_binned, cells_pre_binned_charge_mean(1 + end - baseline_win_binned : end), cells_pre_binned_charge_SEM(1 + end - baseline_win_binned : end), 'color', data_color);
    scatter(post_plot_pts_binned, cells_post_binned_charge_mean, 80, data_color, 'filled');
    errbar(post_plot_pts_binned, cells_post_binned_charge_mean, cells_post_binned_charge_SEM, 'color', data_color);
    rectangle('Position',[depol_range(1) Y_lim_charge*0.47 depol_range(2)-depol_range(1) Y_lim_charge*.03], 'facecolor', [0 0 0]);
    text(depol_range(1) + (depol_range(2)-depol_range(1))/2, Y_lim_charge*.55, 'depol', 'color', [0 0 0], 'horizontalalignment', 'center', 'fontsize', 10, 'fontweight', 'bold');
    text(1000000, 10000, ['n = ', num2str(numRecs)], 'color', data_color);
    xlabel('Time (sec)'); ylabel('Charge transfer (pC)');
    xlim(X_lims); ylim([0 Y_lim_charge]);
    set(gca,'xtick', tick_range_min : tick_range_step : tick_range_max, 'xticklabel', {'-20','-10', '0', '10', '20', '30', '40', '50'})
    
    % Plot the avg normalized binned frequency (across cells)
    
    subplot(2,3,4); hold on;
    title('Normalized Frequency');
    scatter(baseline_plot_pts_binned, norm_cells_pre_binned_freq_mean(1 + end - baseline_win_binned : end), 80, data_color, 'filled');
    errbar(baseline_plot_pts_binned, norm_cells_pre_binned_freq_mean(1 + end - baseline_win_binned : end), norm_cells_pre_binned_freq_SEM(1 + end - baseline_win_binned : end), 'color', data_color);
    scatter(post_plot_pts_binned, norm_cells_post_binned_freq_mean, 80, data_color, 'filled');
    errbar(post_plot_pts_binned, norm_cells_post_binned_freq_mean, norm_cells_post_binned_freq_SEM, 'color', data_color);
    rectangle('Position',[depol_range(1) 125 depol_range(2)-depol_range(1) 4], 'facecolor', [0 0 0]);
    text(depol_range(1) + (depol_range(2)-depol_range(1))/2, 136, 'depol', 'color', [0 0 0], 'horizontalalignment', 'center', 'fontsize', 10, 'fontweight', 'bold');
    text(1000000, 25, ['n = ', num2str(numRecs)], 'color', data_color);
    xlabel('Time (sec)'); ylabel('Frequency (% baseline)');
    xlim(X_lims); ylim([0 150]);
    set(gca,'xtick', tick_range_min : tick_range_step : tick_range_max, 'xticklabel', {'-20','-10', '0', '10', '20', '30', '40', '50'})
    
    % Plot the avg normalized binned amplitude (across cells)
    
    subplot(2,3,5); hold on;
    title('Normalized Amplitude');
    scatter(baseline_plot_pts_binned, norm_cells_pre_binned_amp_mean(1 + end - baseline_win_binned : end), 80, data_color, 'filled');
    errbar(baseline_plot_pts_binned, norm_cells_pre_binned_amp_mean(1 + end - baseline_win_binned : end), norm_cells_pre_binned_amp_SEM(1 + end - baseline_win_binned : end), 'color', data_color);
    scatter(post_plot_pts_binned, norm_cells_post_binned_amp_mean, 80, data_color, 'filled');
    errbar(post_plot_pts_binned, norm_cells_post_binned_amp_mean, norm_cells_post_binned_amp_SEM, 'color', data_color);
    rectangle('Position',[depol_range(1) 125 depol_range(2)-depol_range(1) 4], 'facecolor', [0 0 0]);
    text(depol_range(1) + (depol_range(2)-depol_range(1))/2, 136, 'depol', 'color', [0 0 0], 'horizontalalignment', 'center', 'fontsize', 10, 'fontweight', 'bold');
    text(1000000, 25, ['n = ', num2str(numRecs)], 'color', data_color);
    xlabel('Time (sec)'); ylabel('Amplitude (% baseline)');
    xlim(X_lims); ylim([0 150]);
    set(gca,'xtick', tick_range_min : tick_range_step : tick_range_max, 'xticklabel', {'-20','-10', '0', '10', '20', '30', '40', '50'})
    
    % Plot the avg binned charge transfer (across cells)

    subplot(2,3,6); hold on;
    title('Normalized Charge');
    scatter(baseline_plot_pts_binned, norm_cells_pre_binned_charge_mean(1 + end - baseline_win_binned : end), 80, data_color, 'filled');
    errbar(baseline_plot_pts_binned, norm_cells_pre_binned_charge_mean(1 + end - baseline_win_binned : end), norm_cells_pre_binned_charge_SEM(1 + end - baseline_win_binned : end), 'color', data_color);
    scatter(post_plot_pts_binned, norm_cells_post_binned_charge_mean, 80, data_color, 'filled');
    errbar(post_plot_pts_binned, norm_cells_post_binned_charge_mean, norm_cells_post_binned_charge_SEM, 'color', data_color);
    rectangle('Position',[depol_range(1) 125 depol_range(2)-depol_range(1) 4], 'facecolor', [0 0 0]);
    text(depol_range(1) + (depol_range(2)-depol_range(1))/2, 136, 'depol', 'color', [0 0 0], 'horizontalalignment', 'center', 'fontsize', 10, 'fontweight', 'bold');
    text(1000000, 25, ['n = ', num2str(numRecs)], 'color', data_color);
    xlabel('Time (sec)'); ylabel('Charge transfer (% baseline)');
    xlim(X_lims); ylim([0 150]);
    set(gca,'xtick', tick_range_min : tick_range_step : tick_range_max, 'xticklabel', {'-20','-10', '0', '10', '20', '30', '40', '50'})
    
end

%%
% Calculate summary statistics for freq, amp, and charge
% Compared values: pre vs post1 vs post2

% Calculate mean and SEM for pre, post1, post2

freq_pre_end = mean(cells_pre_binned_freq(:, end - bar_bin_size + 1 : end), 2);
freq_post_start = mean(cells_post_binned_freq(:, 1:bar_bin_size), 2);
freq_post_end = mean(cells_post_binned_freq(:, post_end_t/bin_size : post_end_t/bin_size + bar_bin_size - 1), 2);

freq_pre_end_mean = mean(freq_pre_end, 1);
freq_post_start_mean = mean(freq_post_start, 1);
freq_post_end_mean = mean(freq_post_end, 1);

freq_pre_end_SEM = nanstd(freq_pre_end) ./ sqrt(sum(~isnan(freq_pre_end)));
freq_post_start_SEM = nanstd(freq_post_start) ./ sqrt(sum(~isnan(freq_post_start)));
freq_post_end_SEM = nanstd(freq_post_end) ./ sqrt(sum(~isnan(freq_post_end)));

amp_pre_end = mean(cells_pre_binned_amp(:, end - bar_bin_size + 1 : end), 2);
amp_post_start = mean(cells_post_binned_amp(:, 1:bar_bin_size), 2);
amp_post_end = mean(cells_post_binned_amp(:, post_end_t/bin_size : post_end_t/bin_size + bar_bin_size - 1), 2);

amp_pre_end_mean = mean(amp_pre_end, 1);
amp_post_start_mean = mean(amp_post_start, 1);
amp_post_end_mean = mean(amp_post_end, 1);

amp_pre_end_SEM = nanstd(amp_pre_end) ./ sqrt(sum(~isnan(amp_pre_end)));
amp_post_start_SEM = nanstd(amp_post_start) ./ sqrt(sum(~isnan(amp_post_start)));
amp_post_end_SEM = nanstd(amp_post_end) ./ sqrt(sum(~isnan(amp_post_end)));

charge_pre_end = mean(cells_pre_binned_charge(:, end - bar_bin_size + 1 : end), 2);
charge_post_start = mean(cells_post_binned_charge(:, 1:bar_bin_size), 2);
charge_post_end = mean(cells_post_binned_charge(:, post_end_t/bin_size : post_end_t/bin_size + bar_bin_size - 1), 2);

charge_pre_end_mean = mean(charge_pre_end, 1);
charge_post_start_mean = mean(charge_post_start, 1);
charge_post_end_mean = mean(charge_post_end, 1);

charge_pre_end_SEM = nanstd(charge_pre_end) ./ sqrt(sum(~isnan(charge_pre_end)));
charge_post_start_SEM = nanstd(charge_post_start) ./ sqrt(sum(~isnan(charge_post_start)));
charge_post_end_SEM = nanstd(charge_post_end) ./ sqrt(sum(~isnan(charge_post_end)));

% Calculate the normalized mean and SEM for pre, post1, post2

norm_freq_pre_end = mean(norm_cells_pre_binned_freq(:, end - bar_bin_size + 1 : end), 2);
norm_freq_post_start = mean(norm_cells_post_binned_freq(:, 1:bar_bin_size), 2);
norm_freq_post_end = mean(norm_cells_post_binned_freq(:, post_end_t/bin_size : post_end_t/bin_size + bar_bin_size - 1), 2);

norm_freq_pre_end_mean = mean(norm_freq_pre_end, 1);
norm_freq_post_start_mean = mean(norm_freq_post_start, 1);
norm_freq_post_end_mean = mean(norm_freq_post_end, 1);

norm_freq_pre_end_SEM = nanstd(norm_freq_pre_end) ./ sqrt(sum(~isnan(norm_freq_pre_end)));
norm_freq_post_start_SEM = nanstd(norm_freq_post_start) ./ sqrt(sum(~isnan(norm_freq_post_start)));
norm_freq_post_end_SEM = nanstd(norm_freq_post_end) ./ sqrt(sum(~isnan(norm_freq_post_end)));

norm_amp_pre_end = mean(norm_cells_pre_binned_amp(:, end - bar_bin_size + 1 : end), 2);
norm_amp_post_start = mean(norm_cells_post_binned_amp(:, 1:bar_bin_size), 2);
norm_amp_post_end = mean(norm_cells_post_binned_amp(:, post_end_t/bin_size : post_end_t/bin_size + bar_bin_size - 1), 2);

norm_amp_pre_end_mean = mean(norm_amp_pre_end, 1);
norm_amp_post_start_mean = mean(norm_amp_post_start, 1);
norm_amp_post_end_mean = mean(norm_amp_post_end, 1);

norm_amp_pre_end_SEM = nanstd(norm_amp_pre_end) ./ sqrt(sum(~isnan(norm_amp_pre_end)));
norm_amp_post_start_SEM = nanstd(norm_amp_post_start) ./ sqrt(sum(~isnan(norm_amp_post_start)));
norm_amp_post_end_SEM = nanstd(norm_amp_post_end) ./ sqrt(sum(~isnan(norm_amp_post_end)));

norm_charge_pre_end = mean(norm_cells_pre_binned_charge(:, end - bar_bin_size + 1 : end), 2);
norm_charge_post_start = mean(norm_cells_post_binned_charge(:, 1:bar_bin_size), 2);
norm_charge_post_end = mean(norm_cells_post_binned_charge(:, post_end_t/bin_size : post_end_t/bin_size + bar_bin_size - 1), 2);

norm_charge_pre_end_mean = mean(norm_charge_pre_end, 1);
norm_charge_post_start_mean = mean(norm_charge_post_start, 1);
norm_charge_post_end_mean = mean(norm_charge_post_end, 1);

norm_charge_pre_end_SEM = nanstd(norm_charge_pre_end) ./ sqrt(sum(~isnan(norm_charge_pre_end)));
norm_charge_post_start_SEM = nanstd(norm_charge_post_start) ./ sqrt(sum(~isnan(norm_charge_post_start)));
norm_charge_post_end_SEM = nanstd(norm_charge_post_end) ./ sqrt(sum(~isnan(norm_charge_post_end)));
    
if PLOT_SUMMARY
    
    % Plot the summary statistics for frequency (pre vs post1 vs post2)
    
    figure('units', 'normalized', 'position', [1.01 .15 .98 .7]); hold on;
    subplot(2,3,1); hold on;

    freq_box_h = boxplot(horzcat(freq_pre_end, freq_post_start, freq_post_end));
    set(freq_box_h(:, 1), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    set(freq_box_h(:, 2), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    set(freq_box_h(:, 3), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    scatter(ones(numRecs, 1), freq_pre_end, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    scatter(2 * ones(numRecs, 1), freq_post_start, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    scatter(3 * ones(numRecs, 1), freq_post_end, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    for i = 1:length(freq_pre_end)
        plot([1, 2], [freq_pre_end(i), freq_post_start(i)], 'color', [0 0 0]);
        plot([2, 3], [freq_post_start(i), freq_post_end(i)], 'color', [0 0 0]);
    end
    errbar(1, freq_pre_end_mean, freq_pre_end_SEM, 'color', data_color);
    errbar(2, freq_post_start_mean, freq_post_start_SEM, 'color', data_color);
    errbar(3, freq_post_end_mean, freq_post_end_SEM, 'color', data_color);
    text(2.8, 1, ['n = ', num2str(numRecs)], 'color', data_color);
    title('Frequency');
    ylim([0 Y_lim_freq]);
    set(gca,'xtick', 1:3,'xticklabel',{'baseline', 'early', 'late'});
    ylabel('Frequency (hz)');
    
    [p, h] = signrank(freq_pre_end, freq_post_start); text(1.25, Y_lim_freq*0.9, ['p = ', sprintf('%.4f', p)], 'color', [0 0 0]);
    [p, h] = signrank(freq_post_start, freq_post_end); text(2.25, Y_lim_freq*0.9, ['p = ', sprintf('%.4f', p)], 'color', [0 0 0]);
    
    % Plot the summary statistics for amplitude (pre vs post1 vs post2)
    subplot(2,3,2); hold on;

    amp_box_h = boxplot(horzcat(amp_pre_end, amp_post_start, amp_post_end));
    set(amp_box_h(:, 1), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    set(amp_box_h(:, 2), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    set(amp_box_h(:, 3), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    scatter(ones(numRecs, 1), amp_pre_end, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    scatter(2 * ones(numRecs, 1), amp_post_start, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    scatter(3 * ones(numRecs, 1), amp_post_end, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    for i = 1:length(amp_pre_end)
        plot([1, 2], [amp_pre_end(i), amp_post_start(i)], 'color', [0 0 0]);
        plot([2, 3], [amp_post_start(i), amp_post_end(i)], 'color', [0 0 0]);
    end
    errbar(1, amp_pre_end_mean, amp_pre_end_SEM, 'color', data_color);
    errbar(2, amp_post_start_mean, amp_post_start_SEM, 'color', data_color);
    errbar(3, amp_post_end_mean, amp_post_end_SEM, 'color', data_color);
    text(2.8, 3, ['n = ', num2str(numRecs)], 'color', data_color);
    title('Amplitude');
    ylim([0 Y_lim_amp]);
    set(gca,'xtick', 1:3,'xticklabel',{'baseline', 'early', 'late'});
    ylabel('Amplitude (pA)');
    
    [p, h] = signrank(amp_pre_end, amp_post_start); text(1.25, Y_lim_amp*0.9, ['p = ', sprintf('%.4f', p)], 'color', [0 0 0]);
    [p, h] = signrank(amp_post_start, amp_post_end); text(2.25, Y_lim_amp*0.9, ['p = ', sprintf('%.4f', p)], 'color', [0 0 0]);
    
    % Plot the summary statistics for charge (pre vs post1 vs post2)
    subplot(2,3,3); hold on;

    charge_box_h = boxplot(horzcat(charge_pre_end, charge_post_start, charge_post_end));
    set(charge_box_h(:, 1), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    set(charge_box_h(:, 2), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    set(charge_box_h(:, 3), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    scatter(ones(numRecs, 1), charge_pre_end, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    scatter(2 * ones(numRecs, 1), charge_post_start, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    scatter(3 * ones(numRecs, 1), charge_post_end, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    for i = 1:length(charge_pre_end)
        plot([1, 2], [charge_pre_end(i), charge_post_start(i)], 'color', [0 0 0]);
        plot([2, 3], [charge_post_start(i), charge_post_end(i)], 'color', [0 0 0]);
    end
    errbar(1, charge_pre_end_mean, charge_pre_end_SEM, 'color', data_color);
    errbar(2, charge_post_start_mean, charge_post_start_SEM, 'color', data_color);
    errbar(3, charge_post_end_mean, charge_post_end_SEM, 'color', data_color);
    text(2.8, .8, ['n = ', num2str(numRecs)], 'color', data_color);
    title('Collapsed Charge');
    ylim([0 Y_lim_charge]);
    set(gca,'xtick', 1:3,'xticklabel',{'baseline', 'early', 'late'});
    ylabel('Collapsed Charge (pC)');
    
    [p, h] = signrank(charge_pre_end, charge_post_start); text(1.25, Y_lim_charge*0.9, ['p = ', sprintf('%.4f', p)], 'color', [0 0 0]);
    [p, h] = signrank(charge_post_start, charge_post_end); text(2.25, Y_lim_charge*0.9, ['p = ', sprintf('%.4f', p)], 'color', [0 0 0]);
        
    % Plot the nromalized summary statistics for frequency (pre vs post1 vs post2)
    subplot(2,3,4); hold on;
    
    norm_freq_box_h = boxplot(horzcat(norm_freq_pre_end, norm_freq_post_start, norm_freq_post_end));
    set(norm_freq_box_h(:, 1), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    set(norm_freq_box_h(:, 2), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    set(norm_freq_box_h(:, 3), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    scatter(ones(numRecs, 1), norm_freq_pre_end, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    scatter(2 * ones(numRecs, 1), norm_freq_post_start, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    scatter(3 * ones(numRecs, 1), norm_freq_post_end, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    for i = 1:length(norm_freq_pre_end)
        plot([1, 2], [norm_freq_pre_end(i), norm_freq_post_start(i)], 'color', [0 0 0]);
        plot([2, 3], [norm_freq_post_start(i), norm_freq_post_end(i)], 'color', [0 0 0]);
    end
    errbar(1, norm_freq_pre_end_mean, norm_freq_pre_end_SEM, 'color', data_color);
    errbar(2, norm_freq_post_start_mean, norm_freq_post_start_SEM, 'color', data_color);
    errbar(3, norm_freq_post_end_mean, norm_freq_post_end_SEM, 'color', data_color);
    text(2.8, 20, ['n = ', num2str(numRecs)], 'color', data_color);
    title('Normalized Frequency');
    ylim([0 Y_lim_norm]);
    set(gca,'xtick', 1:3,'xticklabel',{'baseline', 'early', 'late'});
    ylabel('Frequency (relative to baseline)');
    
    [p, h] = signrank(norm_freq_pre_end, norm_freq_post_start); text(1.25, 180, ['p = ', sprintf('%.4f', p)], 'color', [0 0 0]);
    [p, h] = signrank(norm_freq_post_start, norm_freq_post_end); text(2.25, 180, ['p = ', sprintf('%.4f', p)], 'color', [0 0 0]);
    
    % Plot the nromalized summary statistics for amplitude (pre vs post1 vs post2)
    subplot(2,3,5); hold on;
    
    norm_amp_box_h = boxplot(horzcat(norm_amp_pre_end, norm_amp_post_start, norm_amp_post_end));
    set(norm_amp_box_h(:, 1), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    set(norm_amp_box_h(:, 2), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    set(norm_amp_box_h(:, 3), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    scatter(ones(numRecs, 1), norm_amp_pre_end, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    scatter(2 * ones(numRecs, 1), norm_amp_post_start, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    scatter(3 * ones(numRecs, 1), norm_amp_post_end, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    for i = 1:length(norm_amp_pre_end)
        plot([1, 2], [norm_amp_pre_end(i), norm_amp_post_start(i)], 'color', [0 0 0]);
        plot([2, 3], [norm_amp_post_start(i), norm_amp_post_end(i)], 'color', [0 0 0]);
    end
    errbar(1, norm_amp_pre_end_mean, norm_amp_pre_end_SEM, 'color', data_color);
    errbar(2, norm_amp_post_start_mean, norm_amp_post_start_SEM, 'color', data_color);
    errbar(3, norm_amp_post_end_mean, norm_amp_post_end_SEM, 'color', data_color);
    text(2.8, 20, ['n = ', num2str(numRecs)], 'color', data_color);
    title('Normalized Amplitude');
    ylim([0 Y_lim_norm]);
    set(gca,'xtick', 1:3,'xticklabel',{'baseline', 'early', 'late'});
    ylabel('Amplitude (relative to baseline)');
    
    [p, h] = signrank(norm_amp_pre_end, norm_amp_post_start); text(1.25, 180, ['p = ', sprintf('%.4f', p)], 'color', [0 0 0]);
    [p, h] = signrank(norm_amp_post_start, norm_amp_post_end); text(2.25, 180, ['p = ', sprintf('%.4f', p)], 'color', [0 0 0]);
    
    % Plot the nromalized summary statistics for charge transfer (pre vs post1 vs post2)
    subplot(2,3,6); hold on;
 
    norm_charge_box_h = boxplot(horzcat(norm_charge_pre_end, norm_charge_post_start, norm_charge_post_end));
    set(norm_charge_box_h(:, 1), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    set(norm_charge_box_h(:, 2), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    set(norm_charge_box_h(:, 3), 'color', data_color, 'linewidth', 2, 'linestyle', '-');
    scatter(ones(numRecs, 1), norm_charge_pre_end, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    scatter(2 * ones(numRecs, 1), norm_charge_post_start, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    scatter(3 * ones(numRecs, 1), norm_charge_post_end, [], data_color, 'sizedata', 10^2, 'linewidth', 2);
    for i = 1:length(norm_charge_pre_end)
        plot([1, 2], [norm_charge_pre_end(i), norm_charge_post_start(i)], 'color', [0 0 0]);
        plot([2, 3], [norm_charge_post_start(i), norm_charge_post_end(i)], 'color', [0 0 0]);
    end
    errbar(1, norm_charge_pre_end_mean, norm_charge_pre_end_SEM, 'color', data_color);
    errbar(2, norm_charge_post_start_mean, norm_charge_post_start_SEM, 'color', data_color);
    errbar(3, norm_charge_post_end_mean, norm_charge_post_end_SEM, 'color', data_color);
    text(2.8, 20, ['n = ', num2str(numRecs)], 'color', data_color);
    title('Normalized Collapsed Charge');
    ylim([0 Y_lim_norm]);
    set(gca,'xtick', 1:3,'xticklabel',{'baseline', 'early', 'late'});
    ylabel('Collapsed Charge (relative to baseline)');
    
    [p, h] = signrank(norm_charge_pre_end, norm_charge_post_start); text(1.25, 180, ['p = ', sprintf('%.4f', p)], 'color', [0 0 0]);
    [p, h] = signrank(norm_charge_post_start, norm_charge_post_end); text(2.25, 180, ['p = ', sprintf('%.4f', p)], 'color', [0 0 0]);
    
end

%% Save the main data structures into a superstructure

all_data.pre_freq = cells_pre_binned_freq;
all_data.pre_amp = cells_pre_binned_amp;
all_data.pre_charge = cells_pre_binned_charge;
all_data.pre_freq_mean = cells_pre_binned_freq_mean;
all_data.pre_amp_mean = cells_pre_binned_amp_mean;
all_data.pre_charge_mean = cells_pre_binned_charge_mean;
all_data.pre_freq_SEM = cells_pre_binned_freq_SEM;
all_data.pre_amp_SEM = cells_pre_binned_amp_SEM;
all_data.pre_charge_SEM = cells_pre_binned_charge_SEM;

all_data.pre_norm_freq = norm_cells_pre_binned_freq;
all_data.pre_norm_amp = norm_cells_pre_binned_amp;
all_data.pre_norm_charge = norm_cells_pre_binned_charge;
all_data.pre_norm_freq_mean = norm_cells_pre_binned_freq_mean;
all_data.pre_norm_amp_mean = norm_cells_pre_binned_amp_mean;
all_data.pre_norm_charge_mean = norm_cells_pre_binned_charge_mean;
all_data.pre_norm_freq_SEM = norm_cells_pre_binned_freq_SEM;
all_data.pre_norm_amp_SEM = norm_cells_pre_binned_amp_SEM;
all_data.pre_norm_charge_SEM = norm_cells_pre_binned_charge_SEM;

all_data.post_freq = cells_post_binned_freq;
all_data.post_amp = cells_post_binned_amp;
all_data.post_charge = cells_post_binned_charge;
all_data.post_freq_mean = cells_post_binned_freq_mean;
all_data.post_amp_mean = cells_post_binned_amp_mean;
all_data.post_charge_mean = cells_post_binned_charge_mean;
all_data.post_freq_SEM = cells_post_binned_freq_SEM;
all_data.post_amp_SEM = cells_post_binned_amp_SEM;
all_data.post_charge_SEM = cells_post_binned_charge_SEM;

all_data.post_norm_freq = norm_cells_post_binned_freq;
all_data.post_norm_amp = norm_cells_post_binned_amp;
all_data.post_norm_charge = norm_cells_post_binned_charge;
all_data.post_norm_freq_mean = norm_cells_post_binned_freq_mean;
all_data.post_norm_amp_mean = norm_cells_post_binned_amp_mean;
all_data.post_norm_charge_mean = norm_cells_post_binned_charge_mean;
all_data.post_norm_freq_SEM = norm_cells_post_binned_freq_SEM;
all_data.post_norm_amp_SEM = norm_cells_post_binned_amp_SEM;
all_data.post_norm_charge_SEM = norm_cells_post_binned_charge_SEM;
