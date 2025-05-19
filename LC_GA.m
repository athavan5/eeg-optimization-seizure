%[seizure_start_times, seizure_end_times] = LC_predict_seizure_v5();

% Call the function for the first file (e.g., chb01_01.edf)
seizures = predict_seizure(18);

% Access results for a specific channel
disp("predict seizures: ")
disp(seizures)
channel_index = 5; % For channel 5
start_times = seizures{channel_index, 1};
end_times = seizures{channel_index, 2};

disp("Start Times:");
disp(start_times);
disp("End Times:");
disp(end_times);

ground_truth = load('411_Seizure_Times.txt');

% Initial parameter guesses
initial_params = [5, 2, 2.1, 0.5, 5]; % Initial guesses for [window_length, step_size, stdev, drop_ratio, post_peak_window]

% Lower and upper bounds for parameters
lb = [1, 1, 0.5, 0.1, 1]; % Minimum values
ub = [10, 5, 5, 0.9, 10]; % Maximum values

% Optimization
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
optimized_params = fmincon(@(params) evaluate_parameters(params, file_number, ground_truth), ...
                           initial_params, [], [], [], [], lb, ub, [], options);

function loss = evaluate_parameters(params, file_number, ground_truth)
    % Parse parameters
    window_length = params(1);
    step_size = params(2);
    stdev = params(3);
    drop_ratio = params(4);
    post_peak_window = params(5);

    % Run prediction with the given parameters
    seizures = predict_seizure(file_number, window_length, step_size, stdev, drop_ratio, post_peak_window);

    % Compare with ground truth
    % Compute metrics like false positives, false negatives, and total loss
    loss = compute_loss(seizures, ground_truth);
end

function loss = objective_function(seizures, ground_truth)
    %true positive = seizure predictions overlapping
    %false positive = seizure not overlapping
    %false negative = seizures missed
    loss = falsepos+falseneg;
end 