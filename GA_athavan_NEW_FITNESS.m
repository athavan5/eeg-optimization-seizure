% Load and preprocess EEG data
[tt, info, seizures] = load_and_preprocess_eeg('../../Downloads/chb03_01.edf');
channel_names = cellstr(info.SignalLabels);

% Parameters
n = 10; % Population size
c = size(tt, 2); % Number of channels
split = floor(c/2); % Dynamic split point
cc = 0.70; % Crossover chance
mc = 0.005; % Mutation chance
num_elites = 2;
num_gens = 10; % Number of generations

% Create initial population
X = round(rand(c, n));
new_X = zeros(size(X));

% Main GA loop
for gen = 1:num_gens
    % Calculate fitness scores and accuracies
    fitness_scores = zeros(n, 1);
    accuracies = zeros(n, 1);
    for i = 1:n
        [fitness_scores(i), accuracies(i)] = fitness(X(:, i), seizures);
    end
    
    % Non-dominated sorting
    [fronts, crowding_distances] = non_dominated_sort(fitness_scores);
    
    % Selection and elitism
    [new_X, new_fitness_scores, new_accuracies] = selection(X, fitness_scores, accuracies, fronts, crowding_distances, num_elites);
    
    % Crossover
    new_X = crossover(new_X, cc, split, num_elites);
    
    % Mutation
    new_X = mutation(new_X, mc, num_elites);
    
    % Update population
    X = new_X;
    
    % Display progress
    disp(['Generation ', num2str(gen), ' completed']);
end

% Display final results
[fronts, ~] = non_dominated_sort(fitness_scores);
pareto_front = X(:, fronts{1});
pareto_fitness = fitness_scores(fronts{1});
pareto_accuracies = accuracies(fronts{1});

disp('Pareto Front Solutions:');
for i = 1:size(pareto_front, 2)
    solution = pareto_front(:, i);
    selected_channels = find(solution);
    disp(['Solution ', num2str(i), ':']);
    disp(channel_names(selected_channels));
    disp('---');
end

disp('Corresponding Fitness Scores:');
disp(pareto_fitness);

disp('Corresponding Accuracies:');
disp(pareto_accuracies);

% Helper functions
function [tt, info, seizures] = load_and_preprocess_eeg(file_path)
    tt = edfread(file_path);
    info = edfinfo(file_path);
    
    % Apply bandpass filter
    Fs = info.NumSamples(1); % Assuming uniform sampling rate
    [b, a] = butter(4, [0.5 50] / (Fs / 2), 'bandpass');
    tt{:, :} = cellfun(@(x) filtfilt(b, a, x), tt{:, :}, 'UniformOutput', false);
    
    % Get seizure predictions
    seizures = predict_seizure(1); % Assuming file number 1
end

function [fronts, crowding_distances] = non_dominated_sort(fitness_scores)
    % Implement non-dominated sorting algorithm
    % This is a simplified version - consider using a more efficient implementation

    n = size(fitness_scores, 1);
    fronts = {[]};
    dominated = cell(n, 1);
    domination_count = zeros(n, 1);
    
    for i = 1:n
        for j = i+1:n
            if dominates(fitness_scores(i, :), fitness_scores(j, :))
                dominated{i} = [dominated{i}, j];
                domination_count(j) = domination_count(j) + 1;
            elseif dominates(fitness_scores(j, :), fitness_scores(i, :))
                dominated{j} = [dominated{j}, i];
                domination_count(i) = domination_count(i) + 1;
            end
        end
        
        if domination_count(i) == 0
            fronts{1} = [fronts{1}, i];
        end
    end
    
    k = 1;
    while ~isempty(fronts{k})
        Q = [];
        for i = fronts{k}
            for j = dominated{i}
                domination_count(j) = domination_count(j) - 1;
                if domination_count(j) == 0
                    Q = [Q, j];
                end
            end
        end
        k = k + 1;
        fronts{k} = Q;
    end
    
    if isempty(fronts{end})
        fronts(end) = [];
    end
    
    crowding_distances = calculate_crowding_distance(fitness_scores, fronts);
end

function accuracy = calculate_accuracy(actual_seizure_times, predicted_start_times, predicted_end_times)
    % Implement your accuracy calculation here
    % This is a placeholder function - replace with your actual accuracy calculation
    total_seizures = length(actual_seizure_times);
    correct_predictions = 0;
    
    for i = 1:total_seizures
        for j = 1:length(predicted_start_times)
            if abs(actual_seizure_times(i) - predicted_start_times(j)) < 30 % 30-second tolerance
                correct_predictions = correct_predictions + 1;
                break;
            end
        end
    end
    
    accuracy = correct_predictions / total_seizures;
end

function [new_population, new_fitness, new_accuracies] = selection(population, fitness_scores, accuracies, fronts, crowding_distances, num_elites)
    % Implement selection based on non-dominated sorting and crowding distance
    n = size(population, 2);
    new_population = zeros(size(population));
    new_fitness = zeros(size(fitness_scores));
    new_accuracies = zeros(size(accuracies));
    
    % Adjust num_elites based on the size of the first front
    num_elites = min(num_elites, length(fronts{1}));
    
    % Elitism
    elite_indices = fronts{1}(1:num_elites);
    new_population(:, 1:num_elites) = population(:, elite_indices);
    new_fitness(1:num_elites) = fitness_scores(elite_indices);
    new_accuracies(1:num_elites) = accuracies(elite_indices);
    
    % Tournament selection for the rest
    for i = num_elites+1:n
        [winner, winner_fitness, winner_accuracy] = tournament_selection(population, fitness_scores, accuracies, fronts, crowding_distances);
        new_population(:, i) = winner;
        new_fitness(i) = winner_fitness;
        new_accuracies(i) = winner_accuracy;
    end
end

function [winner, winner_fitness, winner_accuracy] = tournament_selection(population, fitness_scores, accuracies, fronts, crowding_distances)
    % Implement tournament selection
    tournament_size = 2;
    candidates = randi(length(crowding_distances), 1, tournament_size);
    
    % Ensure candidates are within bounds
    valid_candidates = all(1 <= candidates & candidates <= length(crowding_distances));
    
    if valid_candidates
        [~, best_candidate] = min([crowding_distances(candidates(1)), crowding_distances(candidates(2))]);
        winner = population(:, candidates(best_candidate));
        winner_fitness = fitness_scores(candidates(best_candidate));
        winner_accuracy = accuracies(candidates(best_candidate));
    else
        % Handle invalid candidates (e.g., select a random individual)
        random_index = randi(size(population, 2));
        winner = population(:, random_index);
        winner_fitness = fitness_scores(random_index);
        winner_accuracy = accuracies(random_index);
    end
end


function new_population = crossover(population, cc, split, num_elites)
    % Implement crossover
    n = size(population, 2);
    new_population = population;
    
    for i = num_elites+1:2:n
        if rand < cc
            parent1 = population(:, i);
            parent2 = population(:, i+1);
            
            child1 = [parent1(1:split); parent2(split+1:end)];
            child2 = [parent2(1:split); parent1(split+1:end)];
            
            new_population(:, i) = child1;
            new_population(:, i+1) = child2;
        end
    end
end

function new_population = mutation(population, mc, num_elites)
    % Implement mutation
    [c, n] = size(population);
    for i = num_elites+1:n
        for j = 1:c
            if rand < mc
                population(j, i) = 1 - population(j, i);
            end
        end
    end
    new_population = population;
end

function result = dominates(a, b)
    result = all(a <= b) && any(a < b);
end

function distances = calculate_crowding_distance(fitness_scores, fronts)
    n = size(fitness_scores, 1);
    distances = zeros(1, n);
    
    for k = 1:length(fronts)
        front = fronts{k};
        if isempty(front)
            continue;
        end
        
        front_size = length(front);
        if front_size == 1
            distances(front) = Inf;
            continue;
        end
        
        for obj = 1:size(fitness_scores, 2)
            [sorted_fitness, sorted_indices] = sort(fitness_scores(front, obj));
            distances(front(sorted_indices(1))) = Inf;
            distances(front(sorted_indices(end))) = Inf;
            for i = 2:front_size-1
                if max(sorted_fitness) - min(sorted_fitness) == 0
                    distances(front(sorted_indices(i))) = Inf;
                else
                    distances(front(sorted_indices(i))) = distances(front(sorted_indices(i))) + ...
                        (sorted_fitness(i+1) - sorted_fitness(i-1)) / (max(sorted_fitness) - min(sorted_fitness));
                end
            end
        end
    end
end

function [fitness_score, accuracy] = fitness(x, seizures)
    % Extract actual seizure information
    actual_seizure_freqs = sum(~cellfun(@isempty, seizures(:,1)));
    
    % Initialize actual seizure times
    actual_seizure_times = [];
    
    % Concatenate non-empty seizure start times
    for i = 1:size(seizures, 1)
        if ~isempty(seizures{i, 1})
            actual_seizure_times = [actual_seizure_times; seizures{i, 1}(:)]; % Append start times
        end
    end
    
    % Predict seizures for selected channels
    selected_channels = find(x);
    predicted_seizures = seizures(selected_channels, :);
    
    % Extract predicted seizure information
    ps = ~cellfun(@isempty, predicted_seizures(:, 1));
    
    % Initialize predicted start and end times
    pst = [];
    pet = [];
    
    % Concatenate non-empty predicted seizure start and end times
    for i = 1:length(selected_channels)
        if ps(i)
            pst = [pst; predicted_seizures{i, 1}(:)]; % Append predicted start times
            pet = [pet; predicted_seizures{i, 2}(:)]; % Append predicted end times
        end
    end
    
    % Calculate fitness score using the fitness function from fitness_ethan.m
    fitness_score = fitness_ethan(x, actual_seizure_freqs, actual_seizure_times, ps, pst, pet);
    
    % Calculate accuracy (you may need to implement this based on your specific requirements)
    accuracy = calculate_accuracy(actual_seizure_times, pst, pet);
end

function fitness_score = fitness_ethan(x, seizure_freqs, seizure_times, ps, pst, pet)
    cw = 0.8; % channel weight
    aw = 1 - cw; % accuracy weight
    start_diff = 8; % buffer of alotted time for seizure start time difference (secs)
    duration_diff = 30; % buffer of alotted time for seizure duration time difference (secs)

    x_size = size(x);
    total_channels = x_size(1);

    % Calculate Score based on # of Channels Used (less is better)
    channels_used = sum(x);
    chan_score = channels_used / total_channels;

    % Calculate Score based on Accuracy
    chan_accs = zeros(1, total_channels);

    for i = 1:total_channels
        if x(i) == 1
            % Ensure we do not exceed the bounds of ps
            if i <= length(ps) && ps(i) % Check if there are predicted seizures for this channel
                % Penalty if number of seizures are not predicted
                if length(pst) ~= seizure_freqs
                    chan_accs(i) = 1;
                elseif isempty(pst) && seizure_freqs == 0
                    chan_accs(i) = 0; % Good score if correctly predicts no seizures
                else
                    % Penalty for difference from start time & duration
                    w = 3;
                    p = 1 / channels_used;
                    for s = 1:seizure_freqs
                        if abs(pst(s) - seizure_times(s, 1)) > start_diff
                            chan_accs(i) = chan_accs(i) + p / seizure_freqs / 2;
                        end
                        if abs(pst(s) - pet(s)) > duration_diff
                            chan_accs(i) = chan_accs(i) + p / seizure_freqs / 2;
                        end
                    end
                end
            else
                chan_accs(i) = 1; % If no predicted seizures or out of bounds, apply penalty or adjust logic as needed.
            end
        else
            chan_accs(i) = 0; % No penalty for unused channels.
        end
    end

    acc_score = sum(chan_accs) / channels_used;
    fitness_score = cw * chan_score + aw * acc_score;
end
