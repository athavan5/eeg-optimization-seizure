% Load raw EEG data from .edf file
raw_data = edfread('../../Downloads/chb20_12.edf'); % Replace with the actual file path

% Initialize a new cell array to store scaled data
scaled_data = raw_data; % Copy the timetable to preserve structure

% Scale each 256x1 double array in the timetable
scaling_factor = 10^4;
for col = 1:width(raw_data) % Iterate over all columns (data variables)
    for row = 1:height(raw_data) % Iterate over all rows (time points)
        % Extract the 256x1 double array
        signal_array = raw_data{row, col}{1}; % Access the cell content (256x1 double)
        
        % Apply the scaling factor
        scaled_signal = signal_array * scaling_factor;
        
        % Update the scaled data
        scaled_data{row, col} = {scaled_signal}; % Store back as a cell
    end
end

% Display the scaled timetable (optional)
disp(scaled_data);


info = edfinfo('../../Downloads/chb20_22.edf');

%Define Sampling Parameters
Fs = 256; % Sampling frequency in Hz
window_size = Fs * 8; % 8 seconds * 256 Hz

%% Step 2: Preprocess the Data
% Apply bandpass filter to remove noise (e.g., 0.5-50 Hz)
eeg_data = scaled_data{:, :};
[b, a] = butter(4, [0.5 50] / (Fs / 2), 'bandpass');
filtered_signals = cellfun(@(x) filtfilt(b, a, x), eeg_data, 'UniformOutput', false);

% Segment the data into 8-second windows
window_size = Fs * 8; % 8 seconds
segmented_signals = cell(size(filtered_signals)); % Initialize

for i = 1:length(filtered_signals)
    segmented_signals{i} = buffer(filtered_signals{i}, window_size, 0, 'nodelay');
end

% Flatten segmented signals into a single array
all_data = []; % Concatenated data
for i = 1:length(segmented_signals)
    all_data = [all_data; segmented_signals{i}(:)];
end

%% Step 3: Feature Extraction
% Extract features for each 8-second segment
num_segments = size(all_data, 2);
features = zeros(num_segments, 5); % Example with 5 features
labels = zeros(num_segments, 1); % Placeholder for seizure labels

for i = 1:num_segments
    segment = all_data(:, i);
    
    % Example features
    features(i, 1) = mean(segment); % Mean
    features(i, 2) = std(segment); % Standard deviation
    features(i, 3) = entropy(segment); % Entropy
    features(i, 4) = max(segment) - min(segment); % Range
    features(i, 5) = bandpower(segment, Fs, [1 30]); % Bandpower in 1-30 Hz range
    
    % Assign label (replace this with actual seizure annotations)
    labels(i) = randi([0, 1]); % Random binary labels (0 = non-seizure, 1 = seizure)
end

% Normalize features
features = (features - min(features)) ./ (max(features) - min(features));

%% Step 4: Define the Optimization Problem
% Decision variables: Binary vector for channel selection
num_channels = size(features, 2);

%% Step 5: Apply NSGA-II
% Optimization options
%s = optimoptions('gamultiobj', 'PopulationSize', 100, 'MaxGenerations', 50);
% Optimization options
options = optimoptions('gamultiobj', ...
    'PopulationSize', 100, ...
    'MaxGenerations', 100, ...
    'InitialPopulationMatrix', [ones(1, num_channels); randi([0 1], 9, num_channels)], ...
    'CrossoverFraction', 0.8, ...
    'MutationFcn', @mutationadaptfeasible, ...
    'Display', 'iter');

% Ensure features and labels are properly sized
if size(features, 1) ~= length(labels)
    error('Number of rows in features must match the length of labels.');
end

% Run NSGA-II
[x, fval] = gamultiobj(@(x) objective_function(x, features, labels), ...
                       num_channels, ...
                       [], [], [], [], ...
                       zeros(num_channels, 1), ...
                       ones(num_channels, 1), ...
                       @(x) constraints(x), ...
                       options);

% Display results
disp('Optimization Results:');
disp(fval);

% Objective function for NSGA-II
function f = objective_function(x, features, labels)
    selected_features = features(:, x > 0.5); % Select features for active channels

    % Debugging: Log selected features
    disp(['Selected features: ', mat2str(find(x > 0.5))]);

    % Check if too few features are selected
    if size(selected_features, 2) < 2
        f = [sum(x), 1]; % Penalize, but don't use Inf
        return;
    end
    
    % Train and evaluate classifier
    try
        model = fitcsvm(selected_features, labels); % SVM example
        cv = crossval(model);
        accuracy = 1 - cv.kfoldLoss; % 1 - error rate
    catch
        f = [sum(x), 1]; % Penalize, but don't use Inf
        return;
    end
    
    % Define objectives
    f = [sum(x), 1 - accuracy]; % Minimize number of channels, Maximize classification accuracy
end

% % Objective function for NSGA-II
% function f = objective_function(x, features, labels)
%     selected_features = features(:, x == 1); % Select features for active channels
% 
%     % Debugging: Log selected features
%     disp(['Selected features: ', mat2str(find(x == 1))]);
% 
%     % Check if no features are selected
%     if isempty(selected_features) || size(selected_features, 2) < 1
%         f(1) = Inf; % Penalize for not selecting features
%         f(2) = Inf; % High error for accuracy objective
%         return;
%     end
% 
%     % Train and evaluate classifier
%     model = fitcsvm(selected_features, labels); % SVM example
% 
%     % Ensure enough data for cross-validation
%     if size(selected_features, 1) < 2
%         f(1) = Inf; % Penalize for insufficient data
%         f(2) = Inf; 
%         return;
%     end
% 
%     accuracy = 1 - crossval(model).kfoldLoss; % 1 - error rate
% 
%     % Define objectives
%     f(1) = sum(x); % Minimize number of channels
%     f(2) = 1 - accuracy; % Maximize classification accuracy
% end

% Constraints: Ensure at least one channel is selected
function [c, ceq] = constraints(x)
    c = -sum(x) + 2; % At least one channel
    ceq = [];
end










% % Extract all the EEG data into a single numeric array
% all_data = []; % Initialize an empty array
% for col = 1:width(scaled_data) % Iterate over all columns (data variables)
%     for row = 1:height(scaled_data) % Iterate over all rows (time points)
%         % Extract the 256x1 double array
%         signal_array = scaled_data{row, col}{1}; % Access the 256x1 array
%         all_data = [all_data; signal_array]; % Concatenate into a single array
%     end
% end
% 
% % Segment the Data into 8-Second Epochs
% % % Segment EEG signals into 8-second windows
% segments = buffer(all_data, window_size, 0, 'nodelay');
