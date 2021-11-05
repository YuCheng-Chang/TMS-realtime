%% TODO

% [X] true amplitude, trial by trial data (adapt phastimate_truephase to also output amplitude)
% [X] visualize equivalent filters, optimize filter object selection
% [X] estimate phase error with Zrenner2018 parameters
% [ ] automated artifact rejection, compare results for data with rejected artifacts
% [ ] include graphical output option in phastimate (for debugging and explanation, including as figure)
% [ ] move filter object creation to separate phastimate_create_filter_objects.m
% [ ] refactor Step 4
% [ ] include more subjects (?)

%%

% Known issues/Developtment TODO:
% - peak frequency has a resolution of 0.5Hz, could be improved
% - change create_epochs to create predetermined number of epochs, currently NUMEPOCHS has to match this exactly
% - when creating subplot figures, subplot should refer to current figure to avoid change of current figure by user moving plot targets
% - there is no protection against designing filters that require longer windows than the available epoch length (filter order depends on peak frequency, range is set in PEAK_FREQUENCY_INTERVAL constant)

%% This is the script that generates the data for the figures and demonstrates the code

% Note:
%  - the script uses estimate_SNR.m to fit 1/f noise and determine SNR
%  - circular variance is used as the measure of dispersion (note that circular variance ranges between 0 and 1 and is equal to 1-PLV)
%  - all circular data is in radians
%  - peak frequency is used to individualize filters in phastimate (0.5 Hz resolution)
%  - epoch creation happens when needed and needs to be identical in different parts of the script
%  - Signal Processing and Global Optimization Toolboxes are required
%  - This script has been tested with Matlab 2017b

%% preliminaries

% check for toolboxes
assert(~isempty(which('designfilt')), 'filter design function designfilt.m not found, is the Signal Processing Toolbox installed?')
assert(~isempty(which('range')), 'statistical function range.m not found, is the Statistics and Machine Learning Toolbox installed?')
assert(~isempty(which('ga')), 'genetic algorithm function ga.m not found, is the Global Optimization Toolbox installed?')

% clear variables, close windows, reset paths
clear all; close all; path(pathdef); clc

% switch to current directory and add relative path to phastimate toolbox
cd(fileparts(getfield(matlab.desktop.editor.getActive, 'Filename')))
% addpath('../phastimate_code/')

% circular statistics functions (simplified from circstat toolbox)
ang_mean = @(x) angle(mean(exp(1i*x)));
ang_diff = @(x, y) angle(exp(1i*x)./exp(1i*y));
ang_var = @(x) 1-abs(mean(exp(1i*x)));
%ang_var2dev = @(v) sqrt(2*v); % circstat preferred formula uses angular deviation (bounded from 0 to sqrt(2)) which is sqrt(2*(1-r))
ang_var2dev = @(v) sqrt(-2*log(1-v)); % formula for circular standard deviation is sqrt(-2*ln(r))

%% set constants
allvec = nan(1,10000);
allts = nan(1,10000);
fs = 1000;
eegsample = 1000;
samplerate = 500;
sample = 0;
downsample = eegsample/samplerate;
idx = downsample;
% filter design method for phastimate (order and peak frequency is variable)
design_phastimate_filter = @(ord, freq, fs) designfilt('bandpassfir', 'FilterOrder', ord, 'CutoffFrequency1', freq-1, 'CutoffFrequency2', freq+1, 'SampleRate', fs, 'DesignMethod', 'window');

NUM_EPOCHS = 17;
HILBERTWIN = 128; % this is an appropriate window for alpha at 1000 Hz
PEAK_FREQUENCY_INTERVAL = [8 14];

%% load resting sate data into a master data table 'T'
disp('Loading the library...');
lib = lsl_loadlib();
disp('Resolving an EEG stream...');
result = {}; 

while isempty(result)
    result = lsl_resolve_byprop(lib,'type','EEG');
end

disp('Opening an inlet...');
inlet = lsl_inlet(result{1});

disp('Now receiving data...');
while 1
    [vec,ts] = inlet.pull_sample();
    if downsample == idx
        idx = 1;
        sample = sample + 1;
        allvec(1,sample) = vec(1);
        allts(1,sample) = ts;
%         disp('sample: ')
%         fprintf('%.2f\t\n',sample);
    
    else
        idx = idx + 1;
    end
    
    if sample == 10000
        break
    end
end

data = allvec;
x = size(data);
T = table('RowNames', {'subj1'});
T.data = data;
T.fs = fs * ones(height(T),1);


T = T(1,:);

clear('data');


%% Step 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine alpha peak frequency and signal to noise ratio
% Note:
% - data is not cleaned before this step
% - figure is created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T.peak_frequency = nan(height(T),1);
T.peak_SNR = nan(height(T),1);

figures_snr = {};
figures_snr{1} = figure;

subplot_index = 1;

for row_index = 1:height(T)
    subject = T(row_index,:).Row{1};
    
    epochs = create_epochs_overlapping(T(row_index,:).data, T(row_index,:).fs); % from continuous data
    
    % estimate SNR and plot
    ax = subplot(10,5,subplot_index);
    [peak_frequency, peak_SNR] = estimate_SNR(epochs, T(row_index,:).fs, PEAK_FREQUENCY_INTERVAL, ax);
     
    % save data
    if ~isempty(peak_frequency)
        T(row_index,:).peak_frequency = peak_frequency;
        T(row_index,:).peak_SNR = peak_SNR;
    end
    
    % beautify axes
    legend off
    title(ax, {subject; ax.Title.String}, 'Interpreter', 'none')
    set(ax, 'XTick', [3, 5, 8, 14, 30])
    if (subplot_index < 46 && row_index <= height(T)-5), xlabel(ax,''); end % unless bottom row
    if (mod(subplot_index,5) ~= 1), ylabel(ax,''); end % cunless first column
    
    subplot_index = subplot_index + 1;
    drawnow
    
    if subplot_index > 50 % new page?
        figures_snr{numel(figures_snr)+1} = figure;
        subplot_index = 1;
    end
end

% save figures as PDF
for i = 1:numel(figures_snr)
    set(figures_snr{i},'Renderer','Painters') %export vectorized
    set(figures_snr{i},'PaperPositionMode','manual','PaperType','a2','PaperOrientation','portrait','PaperUnits','normalized','PaperPosition',[0 0 1 1])
    print(figures_snr{i},['figure_snr (page ' num2str(i) ')'], '-dpdf', '-r0')
end

% exlude subjects without a peak or with negative peak SNR

fprintf('\nRemoving %i/%i entries where no peak could be found [ %s]... ', sum(isnan(T.peak_frequency)), height(T), sprintf('%s ', string(T.Row(isnan(T.peak_frequency)))))
T(isnan(T.peak_frequency), :) = [];
fprintf('done')
fprintf('\nRemoving %i/%i entries with a negative SNR [ %s]... ', sum(T.peak_SNR < 0), height(T), sprintf('%s ', string(T.Row(T.peak_SNR < 0))))
T(T.peak_SNR < 0, :) = [];
fprintf('done')
fprintf('\nThe number of remaining datasets is %i\n', height(T))

clear('subplot_index', 'row_index', 'subject', 'ax', 'epochs', 'peak_frequency', 'peak_SNR', 'i')


%% Step 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine true phase and amplitude, as well as variance of "true" phase
% - epochs are recreated from the data when needed to save memory space
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add columns for epoch-by-epoch data
T.epochs_truephase_mean = nan(height(T), NUM_EPOCHS);
T.epochs_truephase_ang_var = nan(height(T), NUM_EPOCHS);
T.epochs_trueamp_mean = nan(height(T), NUM_EPOCHS);
T.epochs_trueamp_var = nan(height(T), NUM_EPOCHS);

figures_truephase = {};
figures_truephase{1} = figure;

subplot_index = 1;

for row_index = 1:height(T)
    subject = T(row_index,:).Row{1};
    
    epochs = create_epochs_overlapping(T(row_index,:).data, T(row_index,:).fs); % from continuous data

    peak_frequency = T(row_index,:).peak_frequency;
    
    % set-up equivalent filter objects for given peak frequency
    filter_objects = {};
    fs = T(row_index,:).fs;
    for ord = [2 3 4 5] % FIR - windowed sinc
        filter_objects = {filter_objects{:} designfilt('bandpassfir', 'FilterOrder', round(ord * (fs/peak_frequency)), 'CutoffFrequency1', peak_frequency-1, 'CutoffFrequency2', peak_frequency+1, 'SampleRate', fs, 'DesignMethod', 'window')};
    end
    for ord = [3 4 5] % FIR - least squares (equiripple is similar)
        filter_objects = {filter_objects{:} designfilt('bandpassfir', 'FilterOrder', round(ord * (fs/peak_frequency)), 'StopbandFrequency1', peak_frequency-4, 'PassbandFrequency1', peak_frequency-1, 'PassbandFrequency2', peak_frequency+1, 'StopbandFrequency2', peak_frequency+4, 'SampleRate', fs, 'DesignMethod', 'ls')};
    end
    for ord = [4 8 12] % IIR - butterworth
        filter_objects = {filter_objects{:} designfilt('bandpassiir', 'FilterOrder', ord, 'HalfPowerFrequency1', peak_frequency-1, 'HalfPowerFrequency2', peak_frequency+1, 'SampleRate', fs, 'DesignMethod', 'butter')};
    end
    for ord = [4 6 8] % IIR - chebychev I
        filter_objects = {filter_objects{:} designfilt('bandpassiir', 'FilterOrder', ord, 'PassbandFrequency1', peak_frequency-1.5, 'PassbandFrequency2', peak_frequency+1.5, 'PassbandRipple', 0.5, 'SampleRate', fs, 'DesignMethod', 'cheby1')};
    end
    for attenuation = [10 20] % IIR - elliptic
        filter_objects = {filter_objects{:} designfilt('bandpassiir', 'StopbandFrequency1', peak_frequency-2, 'PassbandFrequency1', peak_frequency-1, 'PassbandFrequency2', peak_frequency+1, 'StopbandFrequency2', peak_frequency+2, 'StopbandAttenuation1', attenuation, 'PassbandRipple', 0.5, 'StopbandAttenuation2', attenuation, 'SampleRate', fs, 'DesignMethod', 'ellip', 'MatchExactly', 'passband')};
    end    
    
    [truephase_mean, truephase_variance, trueamp_mean, trueamp_variance] = phastimate_truephase(epochs, filter_objects);
    disp(size(truephase_mean));
    T(row_index,:).epochs_truephase_mean = truephase_mean;
    T(row_index,:).epochs_truephase_ang_var = truephase_variance;
    
    T(row_index,:).epochs_trueamp_mean = trueamp_mean;
    T(row_index,:).epochs_trueamp_var = trueamp_variance;
    
    ax = subplot(10,5,subplot_index); hold on
    histogram(ax, rad2deg(ang_var2dev(T(row_index,:).epochs_truephase_ang_var)), 'BinWidth', 1, 'Normalization', 'cdf', 'DisplayStyle', 'stairs');
    plot(ax, [1 1] .* rad2deg(ang_var2dev(quantile(T(row_index,:).epochs_truephase_ang_var, 0.5))), [0 0.5], 'LineWidth', 2, 'Color', 'red')
    xlim(ax, [0 20]); ylim(ax, [0 1]);
    if (subplot_index > 45 || row_index > height(T)-5), xlabel(ax,'circular deviation (deg)'); end % bottom row
    if (mod(subplot_index,5) == 1), ylabel(ax,'cumulative probability'); end % first column
    title(subject, 'Interpreter', 'none')
    
    subplot_index = subplot_index + 1;
    drawnow
    
    if subplot_index > 50
        figures_truephase{numel(figures_truephase)+1} = figure;
        subplot_index = 1;
    end

end

% save figures as PDF
for i = 1:numel(figures_truephase);
    set(figures_truephase{i},'Renderer','Painters') %export vectorized
    set(figures_truephase{i},'PaperPositionMode','manual','PaperType','a2','PaperOrientation','portrait','PaperUnits','normalized','PaperPosition',[0 0 1 1])
    print(figures_truephase{i},['figure_truephase (page ' num2str(i) ')'], '-dpdf', '-r0')
end

clear('subplot_index', 'row_index', 'subject', 'epochs', 'peak_frequency', 'filter_objects', 'fs', 'ord', 'truephase_mean', 'truephase_variance', 'trueamp_mean', 'trueamp_variance', 'ax', 'h', 'i')


%% Step 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine optimized phastimate parameters and resulting estimate
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%add relevant columns
T.optim_window = nan(height(T), 1);
T.optim_filter_ord = nan(height(T), 1);
T.optim_edge = nan(height(T), 1);
T.optim_ar_ord = nan(height(T), 1);
T.optim_fval = nan(height(T), 1);

fprintf('\nNow running genetic algorithm to find optimized phastimate parameters...')
for row_index = 1:height(T) % iterate through entries
    subject = T(row_index,:).Row{1};
    fprintf('\nProcessing %s ... ', subject);
    assert(T(row_index,:).fs == 1000, 'default bounds for genetic optimization algorithm are set for detecting 8-14 Hz alpha assuming a sample rate of 1kHz');

    epochs = create_epochs_overlapping(T(row_index,:).data, T(row_index,:).fs); % from continuous data

    peak_frequency = T(row_index,:).peak_frequency;

    filter_order_range = 100:250;

    filter_objects_by_order = {}; %the index has to correspond to the order for the genetic algorithm
    for ord = filter_order_range
        filter_objects_by_order{ord} = design_phastimate_filter(ord, peak_frequency, T(row_index,:).fs);
    end
    
    bounds_filter_order = [filter_order_range(1) filter_order_range(end)];
    bounds_window = [400 750];
    bounds_edge = [30 120];
    bounds_ar_order = [5 60];

    % the includemask allows optimizing for a subset of epochs
    % it makes sense to exclude epochs that would also be excluded by the
    % real-time system, e.g. if artifacts are detected so as to not optimize
    % for noisy epochs that wouldn't result in a stimulus anyway
    
    % subselect according to truephase variance
    %includemask = T(row_index,:).epochs_truephase_angdev <= quantile(T(row_index,:).epochs_truephase_angdev, 0.5);

    % subselect according to true amplitude
    includemask = T(row_index,:).epochs_trueamp_mean >= quantile(T(row_index,:).epochs_trueamp_mean, 0.5);
    
    [optimal_parameters, ga_output] = phastimate_optimize(epochs(1:ceil(end/2),includemask), T(row_index,:).epochs_truephase_mean(includemask), filter_objects_by_order, bounds_filter_order, bounds_window, bounds_edge, bounds_ar_order, HILBERTWIN);

    % rerun phastimate with the optimized settings to confirm result
    D = design_phastimate_filter(optimal_parameters.filter_order, peak_frequency, T(row_index,:).fs);
    [estphase, estamp] = phastimate(epochs(((-optimal_parameters.window_length+1):0)+ceil(end/2),:), D, optimal_parameters.edge, optimal_parameters.ar_order, 128);
    [phase] =  phase1(optimal_parameters.window_length, D, optimal_parameters.edge, optimal_parameters.ar_order);
    % sanity check if the angular deviation matches the result of the optimization
    phases_error = ang_diff(T(row_index,:).epochs_truephase_mean, estphase);
    assert(abs(optimal_parameters.fval - ang_var(phases_error(includemask))) < 0.01, 'could not replicate result of optimization, were the same filters used?')
    
    %TODO: save output of ga, including number of generations etc.

    T(row_index,:).optim_window = optimal_parameters.window_length;
    T(row_index,:).optim_filter_ord = optimal_parameters.filter_order;
    T(row_index,:).optim_edge = optimal_parameters.edge;
    T(row_index,:).optim_ar_ord = optimal_parameters.ar_order;
    T(row_index,:).optim_fval = optimal_parameters.fval;
    

end
fprintf('\nDone.\n')

clear('row_index', 'subject', 'epochs', 'peak_frequency', 'filter_order_range', 'ord', 'bounds_ar_order', 'filter_objects_by_order', 'bounds_edge', 'bounds_filter_order', 'bounds_window', 'D', 'includemask', 'estamp', 'estphase', 'optimal_parameters', 'ga_output', 'phases_error')

%% save data

fprintf('\nSaving the results table...')
save('results_table', 'T', '-v7.3')
fprintf('\nDone.\n')

% NOTE: It's possible to load the data and continue the script below
T(ismember(T.Row, {'SCREEN_001', 'SCREEN_017', 'SCREEN2_085'}),:) = [];
%%T.epochs_truephase_ang_var = (T.epochs_truephase_angdev.^2)/2;



