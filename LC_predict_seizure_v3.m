% Load raw EEG data from .edf file
raw_data = edfread("chb_files\chb01_03.edf"); % PICK FILE

params.Fs = 256;
Fs = params.Fs;

% Specify the channel to analyze
channel_index = 18; % PICK CHANNEL
channel_label = raw_data.Properties.VariableNames{channel_index};
channel_data = raw_data.(channel_label);
full_signal = cell2mat(channel_data);
total_samples = length(full_signal);

% Setting windows
window_length = 5; % WINDOW DURATION (s)
step_size = 2; % STEP (s)
window_samples = window_length * Fs;
step_samples = step_size * Fs;

% Preallocating feature arrays
num_windows = floor((total_samples - window_samples) / step_samples) + 1;
rms_values = zeros(num_windows, 1);
std_values = zeros(num_windows, 1);
peak_counts = zeros(num_windows, 1);

% Loop through the signal in windows
for i = 1:num_windows
    % Define start and end indices for the window
    start_idx = (i - 1) * step_samples + 1;
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

rms_values = rms_values(1:num_windows);
std_values = std_values(1:num_windows);
peak_counts = peak_counts(1:num_windows);

% Moving baseline for features
window_num = 14; % PICK # WINDOWS (for baseline calculation)
baseline_rms = zeros(size(rms_values));
baseline_peak_density = zeros(size(peak_counts));
baseline_std = zeros(size(std_values));

% Calculates baseline (excluding current window)
for i = 1:num_windows
    if i > window_num % enough past windows available, calculate baseline with prev window_num windows
        baseline_rms(i) = mean(rms_values(i - window_num:i - 4)); % defines indices of prev window_num windows
        baseline_peak_density(i) = mean(peak_counts(i - window_num:i - 4));
        baseline_std(i) = mean(std_values(i - window_num:i - 4));
    else % first couple windows, calculates baseline with all initial windows
        baseline_rms(i) = mean(rms_values(1:i - 1));
        baseline_peak_density(i) = mean(peak_counts(1:i - 1));
        baseline_std(i) = mean(std_values(1:i - 1));
    end
end

% Detect anomalies
stdev = 3; % PICK STANDARD DEVIATION (baseline comparison)
anomalies_rms = rms_values > baseline_rms + stdev * movstd(rms_values, window_num); 
anomalies_peak_density = peak_counts > baseline_peak_density + stdev * movstd(peak_counts, window_num);
anomalies_std = std_values > baseline_std + stdev * movstd(std_values, window_num);

% Convert window indices to time in seconds
time_seconds = (0:num_windows - 1) * step_size;

% Visualization of Baseline Verification
figure;
plot(time_seconds, rms_values, 'b', 'DisplayName', 'RMS');
hold on;
plot(time_seconds, baseline_rms, 'r--', 'DisplayName', 'Baseline RMS');
plot(time_seconds, baseline_rms + stdev * movstd(rms_values, window_num), 'g--', 'DisplayName', 'Upper Threshold');
plot(time_seconds(anomalies_rms), rms_values(anomalies_rms), 'ko', 'DisplayName', 'Anomalies'); % Plot anomaly dots
xlabel('Time (s)');
ylabel('RMS Amplitude');
title('RMS Amplitude with Anomalies');
legend;

% Plot features for analysis
figure;
sgtitle([channel_label, ' features'], 'FontSize', 14, 'FontWeight', 'bold');

% RMS plot
subplot(3, 1, 1);
plot(time_seconds, rms_values, 'b', 'DisplayName', 'RMS');
hold on;
plot(time_seconds, baseline_rms, 'r', 'DisplayName', 'Baseline RMS');
plot(time_seconds, baseline_rms + stdev * movstd(rms_values, window_num), 'r--', 'DisplayName', 'Upper Threshold');
plot(time_seconds(anomalies_rms), rms_values(anomalies_rms), 'ko', 'DisplayName', 'Anomalies'); % Plot anomaly dots
xlabel('Time (s)');
ylabel('RMS Amplitude');
title('RMS Amplitude with Anomalies');
legend;

% Peak density plot
subplot(3, 1, 2);
plot(time_seconds, peak_counts, 'g', 'DisplayName', 'Peak Density');
hold on;
plot(time_seconds, baseline_peak_density, 'r--', 'DisplayName', 'Baseline Peak Density');
plot(time_seconds, baseline_peak_density + stdev * movstd(peak_counts, window_num), 'r--', 'DisplayName', 'Upper Threshold');
plot(time_seconds(anomalies_peak_density), peak_counts(anomalies_peak_density), 'ko', 'DisplayName', 'Anomalies'); % Plot anomaly dots
xlabel('Time (s)');
ylabel('Peak Density');
title('Peak Density with Anomalies');
legend;

% Standard deviation plot
subplot(3, 1, 3);
plot(time_seconds, std_values, 'm', 'DisplayName', 'Standard Deviation');
hold on;
plot(time_seconds, baseline_std, 'r--', 'DisplayName', 'Baseline Std Dev');
plot(time_seconds, baseline_std + stdev * movstd(std_values, window_num), 'r--', 'DisplayName', 'Upper Threshold');
plot(time_seconds(anomalies_std), std_values(anomalies_std), 'ko', 'DisplayName', 'Anomalies'); % Plot anomaly dots
xlabel('Time (s)');
ylabel('Standard Deviation');
title('Standard Deviation with Anomalies');
legend;
