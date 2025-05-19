%%%
% Version 1
% - checks to verify data for each channel
% - plots signal (can choose channel)
% - plots rms values and peak density with time windows

% Load raw EEG data from .edf file
raw_data = edfread("chb_files\chb01_03.edf"); % PICK FILE

params.Fs = 256;
Fs = params.Fs;

%{
%print out data to check channels
% Display the first few rows for each channel
disp('First few rows for different channels:');
for idx = 5:10 % CHANNELS TO CHECK
    channel_label = raw_data.Properties.VariableNames{idx};
    channel_data = raw_data.(channel_label); % Extract first second of data
    channel_data = cell2mat(channel_data);
    fprintf('Channel %s:\n', channel_label);
    disp(channel_data(5:8)); % DISPLAY TIME WINDOW
end
%}

% Specify the channel to analyze
channel_index = 15; % PICK CHANNEL
channel_label = raw_data.Properties.VariableNames{channel_index};
channel_data = raw_data.(channel_label);
full_signal = cell2mat(channel_data);
total_samples = length(full_signal);

%plotting the signal
plotSignal(2900, 3200, full_signal, channel_label, params)

% setting windows
window_length = 8; % WINDOW DURATION (s)
window_samples = window_length * Fs;

% preallocating feature arrays
num_windows = floor(total_samples/window_samples);
rms_values = zeros(num_windows, 1);
peak_counts = zeros(num_windows, 1);
std_values = zeros(num_windows, 1);

% loop through the signal in windows
for i = 1:num_windows
    % Define start and end indices for the window
    start_idx = (i - 1) * window_samples + 1;
    end_idx = i * window_samples;

    % Extract the windowed segment
    segment = full_signal(start_idx:end_idx);

    % Compute features
    rms_values(i) = rms(segment); % Root mean square
    std_values(i) = std(segment); % Variability
    [~, peak_locs] = findpeaks(segment); % Peaks
    peak_counts(i) = length(peak_locs); % Peak density
end

time_windows = linspace(0, total_samples / Fs, num_windows); % Time vector for windows

% Plot features for analysis
figure;
sgtitle([channel_label, ' features'], 'FontSize', 14, 'FontWeight', 'bold');
subplot(3, 1, 1);
plot(time_windows, rms_values);
xlabel('Time (s)');
ylabel('RMS Amplitude');
title('RMS Amplitude Over Time');

subplot(3, 1, 2);
plot(time_windows, peak_counts);
xlabel('Time (s)');
ylabel('Peak Count');
title('Peak Density Over Time');

subplot(3, 1, 3);
plot(time_windows, std_values);
xlabel('Time (s)');
ylabel('Standard Deviation');
title('Standard Deviation Over Time');

%{
% Example: Threshold to Identify Preictal Activity
rms_threshold = 50; % Adjust based on your data
peak_threshold = 20; % Adjust based on your data

preictal_idx = find(rms_values > rms_threshold & peak_counts > peak_threshold, 1);
if ~isempty(preictal_idx)
    disp(['Preictal stage starts at window: ', num2str(preictal_idx), ...
          ', corresponding to time: ', num2str(time_windows(preictal_idx)), ' seconds.']);
else
    disp('No preictal stage detected in this segment.');
end
%}


function plotSignal(start_time, end_time, full_signal, channel_label, params)
    Fs = params.Fs;
    start_sample = round(start_time*Fs);
    end_sample = round(end_time*Fs);

    total_samples = length(full_signal); % Total number of samples
    time_vector = linspace(start_time, end_time, end_sample - start_sample + 1); % Time vector
    signal_to_plot = full_signal(start_sample:end_sample); % Extract corresponding signal segment
    
    % Plot the first second of data for Channel
    figure;
    plot(time_vector, signal_to_plot);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(['EEG Signal - ', channel_label]);
end




% Scale the EEG signals to standardize intensity range
%scaling_factor = 10^4;
%standardized_signal = raw_data * scaling_factor;

%Define Sampling Parameters
%Fs = 256; % Sampling frequency in Hz, records 256 amplitude (voltage) values per second
%window_size = Fs * 8; % 8 seconds * 256 Hz

% Segment EEG signals into 8-second windows
%segments = buffer(standardized_signal, window_size, 0, 'nodelay');
