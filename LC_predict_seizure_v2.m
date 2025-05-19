% Load raw EEG data from .edf file
raw_data = edfread("chb_files\chb01_03.edf"); % PICK FILE

params.Fs = 256;
Fs = params.Fs;

% Specify the channel to analyze
channel_index = 19; % PICK CHANNEL
channel_label = raw_data.Properties.VariableNames{channel_index};
channel_data = raw_data.(channel_label);
full_signal = cell2mat(channel_data);
total_samples = length(full_signal);

%plotting the signal
%plotSignal(2994, 3010, full_signal, channel_label, params)

% setting windows
window_length = 8; % WINDOW DURATION (s)
step_size = 2; % STEP (s)
window_samples = window_length * Fs;
step_samples = step_size*Fs;

% preallocating feature arrays
%num_windows = floor(total_samples/window_samples);
num_windows = floor((total_samples - window_samples) / step_samples) + 1;
rms_values = zeros(num_windows, 1);
std_values = zeros(num_windows, 1);
peak_counts = zeros(num_windows, 1);

% loop through the signal in windows
for i = 1:num_windows
    % Define start and end indices for the window
    start_idx = (i - 1) * window_samples + 1;
    end_idx = start_idx + window_samples - 1;

    % Check if end_idx exceeds total_samples
    if end_idx > total_samples
        break; % Exit the loop if the window goes beyond the signal length
    end

    % Extract the windowed segment
    segment = full_signal(start_idx:end_idx);

    % Compute features
    rms_values(i) = rms(segment); % Root mean square
    std_values(i) = std(segment); % Variability
    [~, peak_locs] = findpeaks(segment); % Peaks
    peak_counts(i) = length(peak_locs); % Peak density
end

rms_values = rms_values(1:i-1);
std_values = std_values(1:i-1);
peak_counts = peak_counts(1:i-1);

% Moving baseline for RMS
baseline_rms = movmean(rms_values, 10); % 5-window rolling mean
baseline_peak_density = movmean(peak_counts, 10); % 5-window rolling mean
baseline_std = movmean(std_values, 10); % 5-window rolling mean

% detect anomalies
anomalies_rms = rms_values > baseline_rms + 2 * movstd(rms_values, 5); % 2 std dev above baseline
anomalies_peak_density = peak_counts > baseline_peak_density + 2 * movstd(peak_counts, 5);
anomalies_std = std_values > baseline_std + 2 * movstd(std_values, 5);

% Adjust time windows or indices to match computed feature array length
actual_num_windows = length(rms_values);
window_indices = 1:actual_num_windows; % X-axis indices for computed windows

% Plot features for analysis
%time_windows = linspace(0, total_samples / Fs, num_windows); % Time vector for windows
figure;
sgtitle([channel_label, ' features'], 'FontSize', 14, 'FontWeight', 'bold');

% RMS plot
subplot(3, 1, 1);
plot(window_indices, rms_values, 'b', 'DisplayName', 'RMS');
hold on;
plot(window_indices, baseline_rms(1:actual_num_windows), 'r', 'DisplayName', 'Baseline RMS');
plot(find(anomalies_rms(1:actual_num_windows)), rms_values(anomalies_rms(1:actual_num_windows)), 'ko', 'DisplayName', 'Anomalies');
xlabel('Window Index');
ylabel('RMS Amplitude');
title('RMS Amplitude with Anomalies');
legend;

% Peak density plot
subplot(3, 1, 2);
plot(window_indices, peak_counts, 'g', 'DisplayName', 'Peak Density');
hold on;
plot(window_indices, baseline_peak_density(1:actual_num_windows), 'r', 'DisplayName', 'Baseline Peak Density');
plot(find(anomalies_peak_density(1:actual_num_windows)), peak_counts(anomalies_peak_density(1:actual_num_windows)), 'ko', 'DisplayName', 'Anomalies');
xlabel('Window Index');
ylabel('Peak Density');
title('Peak Density with Anomalies');
legend;

% Standard deviation plot
subplot(3, 1, 3);
plot(window_indices, std_values, 'm', 'DisplayName', 'Standard Deviation');
hold on;
plot(window_indices, baseline_std(1:actual_num_windows), 'r', 'DisplayName', 'Baseline Std Dev');
plot(find(anomalies_std(1:actual_num_windows)), std_values(anomalies_std(1:actual_num_windows)), 'ko', 'DisplayName', 'Anomalies');
xlabel('Window Index');
ylabel('Standard Deviation');
title('Standard Deviation with Anomalies');
legend;
