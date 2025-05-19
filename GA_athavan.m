% Load and preprocess EEG data
[tt, info] = load_and_preprocess_eeg('chb01_01.edf');

channel_names = {
    'FP1-F7', 'F7-T7', 'T7-P7', 'P7-O1', 'FP1-F3', 'F3-C3', ...
    'C3-P3', 'P3-O1', 'FP2-F4', 'F4-C4', 'C4-P4', 'P4-O2', ...
    'FP2-F8', 'F8-T8', 'T8-P8', 'P8-O2', 'FZ-CZ', 'CZ-PZ', ...
    'P7-T7', 'T7-FT9', 'FT9-FT10', 'FT10-T8', 'T8-P8'
};

% Parameters
n = 50; % Increased population size
c = size(tt, 2); % Number of channels
split = floor(c/2); % Dynamic split point
cc = 0.80; % Crossover chance
mc = 0.005; % Slightly increased mutation chance
num_elites = 2;
num_gens = 3; % Increased number of generations

% Create initial population
X = round(rand(c, n));
new_X = zeros(size(X));

% Main GA loop
for gen = 1:num_gens
    % Calculate fitness scores
    fitness_scores = zeros(n, 2); % Two objectives: channel count and accuracy
    for i = 1:n
        fitness_scores(i, :) = fitness(X(:, i), tt, info);
    end
    
    % Non-dominated sorting
    [fronts, crowding_distances] = non_dominated_sort(fitness_scores);
    
    % Selection and elitism
    [new_X, new_fitness_scores] = selection(X, fitness_scores, fronts, crowding_distances, num_elites);
    
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
pareto_fitness = fitness_scores(fronts{1}, :);
disp('Pareto Front Solutions:');
%disp(pareto_front);
for i = 1:size(pareto_front, 1)
    solution = pareto_front(i, :);
    selected_channels = [];
    for j = 1:length(solution)
        if solution(j) == 1
            start_idx = (j-1)*6 + 1;
            end_idx = min(j*6, 23);
            selected_channels = [selected_channels, start_idx:end_idx];
        end
    end
    
    disp(['Solution ', num2str(i), ':']);
    disp(channel_names(selected_channels));
    disp('---');
end
disp('Corresponding Fitness Scores:');
disp(pareto_fitness);

function [tt, info] = load_and_preprocess_eeg(file_path)
    tt = edfread(file_path);
    info = edfinfo(file_path);
    
    % Apply bandpass filter
    Fs = info.NumSamples(1); % Assuming uniform sampling rate
    [b, a] = butter(4, [0.5 50] / (Fs / 2), 'bandpass');
    tt{:, :} = cellfun(@(x) filtfilt(b, a, x), tt{:, :}, 'UniformOutput', false);
end

function [score] = fitness(x, tt, info)
    selected_channels = find(x);
    num_channels = sum(x);
    
    if num_channels < 2
        score = [1, 1]; % Penalize solutions with too few channels
        return;
    end
    
    % Extract features from selected channels
    features = extract_features(tt(:, selected_channels));
    
    % Calculate accuracy (placeholder - replace with actual seizure prediction)
    accuracy = rand(); % Simulated accuracy
    
    % Normalize objectives
    chan_score = num_channels / size(tt, 2);
    acc_score = 1 - accuracy;
    
    score = [chan_score, acc_score];
end

function features = extract_features(data)
    % Extract relevant features from EEG data
    % This is a placeholder - implement actual feature extraction
    features = zeros(height(data), 5);
    for i = 1:height(data)
        segment = cell2mat(data{i, :});
        features(i, 1) = mean(segment(:));
        features(i, 2) = std(segment(:));
        features(i, 3) = entropy(segment(:));
        features(i, 4) = max(segment(:)) - min(segment(:));
        features(i, 5) = bandpower(segment(:), 256, [1 30]);
    end
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

function [new_population, new_fitness] = selection(population, fitness_scores, fronts, crowding_distances, num_elites)
    % Implement selection based on non-dominated sorting and crowding distance
    if length(crowding_distances) ~= size(population, 2)
        error('Length of crowding_distances does not match population size');
    end

    n = size(population, 2);
    new_population = zeros(size(population));
    new_fitness = zeros(size(fitness_scores));
    
    % Elitism
    elite_indices = fronts{1}(1:num_elites);
    new_population(:, 1:num_elites) = population(:, elite_indices);
    new_fitness(1:num_elites, :) = fitness_scores(elite_indices, :);
    
    % Tournament selection for the rest
    for i = num_elites+1:n
        [winner, winner_fitness] = tournament_selection(population, fitness_scores, fronts, crowding_distances);
        new_population(:, i) = winner;
        new_fitness(i, :) = winner_fitness;
    end
end

function [winner, winner_fitness] = tournament_selection(population, fitness_scores, fronts, crowding_distances)
    % Implement tournament selection
    tournament_size = 2;
    candidates = randi(length(crowding_distances), 1, tournament_size);
    
    % Ensure candidates are within bounds
    valid_candidates = all(1 <= candidates & candidates <= length(crowding_distances));
    
    if valid_candidates
        [~, best_candidate] = min([crowding_distances(candidates(1)), crowding_distances(candidates(2))]);
        winner = population(:, candidates(best_candidate));
        winner_fitness = fitness_scores(candidates(best_candidate), :);
    else
        % Handle invalid candidates (e.g., select a random individual)
        random_index = randi(size(population, 2));
        winner = population(:, random_index);
        winner_fitness = fitness_scores(random_index, :);
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