function seizures = predict_seizure(file_number)
    % Predict seizures for all 23 channels in a given EEG file.
    % Inputs:
    %   file_number - The index of the file to load (e.g., 1 for "chb01_01.edf")
    % Outputs:
    %   seizures - A 3D cell array (23 x 2 x ?) containing seizure start and end times for each channel.

    % Generate the filename based on the input file number
    file_prefix = "chb_files\chb01_";
    file_suffix = sprintf('%02d.edf', file_number);
    filename = file_prefix + file_suffix;

    % Load raw EEG data from .edf file
    raw_data = edfread(filename);
    
    % Sampling frequency
    Fs = 256;
    
    % Window Parameters
    window_length = 5; % WINDOW DURATION (s)
    step_size = 2; % STEP (s)
    window_samples = window_length * Fs;
    step_samples = step_size * Fs;

    % Initialize output structure
    num_channels = 23;
    seizures = cell(num_channels, 2);

    % Loop through each channel
    for channel_index = 1:num_channels
        % Extract channel data
        channel_label = raw_data.Properties.VariableNames{channel_index};
        channel_data = raw_data.(channel_label);
        full_signal = cell2mat(channel_data);
        total_samples = length(full_signal);

        % Preallocating feature arrays
        num_windows = floor((total_samples - window_samples) / step_samples) + 1;
        rms_values = zeros(num_windows, 1);
        %std_values = zeros(num_windows, 1);
        %peak_counts = zeros(num_windows, 1);
        
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
            %std_values(i) = std(segment); % Variability
            %[~, peak_locs] = findpeaks(segment); % Peaks
            %peak_counts(i) = length(peak_locs); % Peak density
        end
        
        rms_values = rms_values(1:num_windows);
        %std_values = std_values(1:num_windows);
        %peak_counts = peak_counts(1:num_windows);
        
        % Moving baseline for features (rolling median)
        window_num = 14; % PICK # WINDOWS (for baseline calculation)
        baseline_rms = movmedian(rms_values, window_num);
        %baseline_peak_density = movmedian(peak_counts, window_num);
        %baseline_std = movmedian(std_values, window_num);
        
        % Calculates baseline using only prior windows
        recent_offset = 5; % Exclude the last 5 windows for baseline calculation
        
        for i = 1:num_windows
            if i > window_num + recent_offset
                % Use windows further back, excluding the recent `recent_offset` windows
                baseline_rms(i) = mean(rms_values(i - window_num - recent_offset:i - recent_offset - 1));
                %baseline_peak_density(i) = mean(peak_counts(i - window_num - recent_offset:i - recent_offset - 1));
                %baseline_std(i) = mean(std_values(i - window_num - recent_offset:i - recent_offset - 1));
            elseif i > recent_offset
                % Use all windows available so far, excluding the recent offset
                baseline_rms(i) = mean(rms_values(1:i - recent_offset - 1));
                %baseline_peak_density(i) = mean(peak_counts(1:i - recent_offset - 1));
                %baseline_std(i) = mean(std_values(1:i - recent_offset - 1));
            else
                % For early cases where offset windows are not available
                baseline_rms(i) = mean(rms_values(1:i - 1));
                %baseline_peak_density(i) = mean(peak_counts(1:i - 1));
                %baseline_std(i) = mean(std_values(1:i - 1));
            end
        end
        
        % Define upper thresholds
        stdev = 2.1; % PICK STANDARD DEVIATION (baseline comparison)
        % Use IQR for a robust threshold
        threshold_rms = baseline_rms + stdev * iqr(rms_values);
        %threshold_peak_density = baseline_peak_density + stdev * iqr(peak_counts);
        %threshold_std = baseline_std + stdev * iqr(std_values);
        
        % Detect anomalies
        anomalies_rms = rms_values > threshold_rms;
        %anomalies_peak_density = peak_counts > threshold_peak_density;
        %anomalies_std = std_values > threshold_std;
        
        %%% SEIZURE start and end times
        
        % Parameters for seizure end detection
        min_seizure_duration = 12; % Minimum seizure duration in seconds
        max_seizure_duration = 300; % Maximum seizure duration in seconds
        drop_ratio = 0.5; % Define drop threshold as 50% below the peak
        post_peak_window = 5; % Search window for a "recovery" peak after a drop
        
        % Convert duration to number of windows
        min_seizure_windows = ceil(min_seizure_duration / step_size);
        max_seizure_windows = floor(max_seizure_duration / step_size);
        post_peak_samples = ceil(post_peak_window / step_size);
        
        % Initialize arrays to store start and end times
        seizure_start_times = [];
        seizure_end_times = [];
        
        in_seizure = false; % State variable
        
        % Iterate through the anomaly array to find start and end times
        for i = 1:num_windows
            if anomalies_rms(i) && ~in_seizure
                % Start of seizure
                seizure_start = i;
                in_seizure = true;
            end
        
            if in_seizure
                % Find peak within the seizure
                current_peak = max(rms_values(seizure_start:i));
        
                % Check for drop condition
                if rms_values(i) < drop_ratio * current_peak
                    % Search for a recovery peak within the next `post_peak_window`
                    recovery_peak = false;
                    for j = i+1:min(i+post_peak_samples, num_windows)
                        if rms_values(j) > rms_values(i)
                            recovery_peak = true;
                            break;
                        end
                    end
        
                    if recovery_peak
                        seizure_end = i;
                        seizure_duration = seizure_end - seizure_start;
        
                        % Validate seizure duration
                        if seizure_duration >= min_seizure_windows && ...
                           seizure_duration <= max_seizure_windows
                            % Save seizure times
                            seizure_start_times = [seizure_start_times, seizure_start];
                            seizure_end_times = [seizure_end_times, seizure_end];
                        end
        
                        % Reset state
                        in_seizure = false;
                    end
                end
            end
        end
        
        % Convert window indices to time in seconds
        seizure_start_times = seizure_start_times * step_size;
        seizure_end_times = seizure_end_times * step_size;
        
        % Store results for the channel
        seizures{channel_index, 1} = seizure_start_times;
        seizures{channel_index, 2} = seizure_end_times;
    end
end
