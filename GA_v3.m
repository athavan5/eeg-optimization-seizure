tt = edfread('chb01_01.edf');
info = edfinfo('chb01_01.edf');

% Parameters
n = 10; % # of initial solutions
c = 23; % # of channels
split = 12; % index of split; channels up to and including are 1st string
cc = 0.70; % crossover chance
mc = 0.005; % mutation chance

num_elites = 2;
num_gens = 10;

% Create random c x n matrix of 0s and 1s (X)
X = round(rand(c,n));
new_X = zeros(size(X));

% Create/Calculate fitness score vector
fitness_scores = zeros(n,1);
for i = 1:n
    fitness_scores(i) = fitness(X(:,i),tt,info);
end

% Create empty matrices of histories
history_X = zeros(c,n,num_gens);
history_mpool = zeros(c,n,num_gens);
history_fscore = zeros(n,num_gens);
history_fsums = zeros(1,num_gens);

% Repeat algorithm for specified number of generations
for gen = 1:num_gens

    % Sort X into new_X based on fitness scores
    for i = 1:n
        max_score_index = 0;
        max_score = -inf;
        
        for j = 1:n
            if fitness_scores(j) > max_score
                max_score_index = j;
                max_score = fitness_scores(j);
            end
        end
        
        new_X(:,i) = X(:,max_score_index);
        fitness_scores(max_score_index) = 0;
    end
    
    % Recalculate fitness scores in new_X using fitness function
    for i = 1:n
        fitness_scores(i) = fitness(new_X(:,i),tt,info);
    end

    % Create Normalized Fitness Scores
    fitness_sum = sum(fitness_scores);
    fitness_norms = fitness_scores/fitness_sum;

    % Mating Pool Creation
    mating_pool = zeros(c,n);
    for i = 1:n
        r = rand;
        p = 0;
        
        for j = 1:n
            if r < p + fitness_norms(j)
                mating_pool(:,i) = new_X(:,j);
                break
            else
                p = p + fitness_norms(j);
            end
        end
    end

    % Crossover
    cross_sols_remaining = n - num_elites;

    for i = num_elites+1:2:n
        if cross_sols_remaining < 2
            break
        end

        crossover_decider = rand;
        rand_index1 = randi(n);
        rand_index2 = randi(n);
        parent1 = mating_pool(:,rand_index1);
        parent2 = mating_pool(:,rand_index2);

        if crossover_decider < cc
            for chan = 1:split
                new_X(chan,i) = parent1(chan);
                new_X(chan,i+1) = parent2(chan);
            end
            for chan = split+1:c
                new_X(chan,i) = parent2(chan);
                new_X(chan,i+1) = parent1(chan);
            end
        else
            new_X(:,i) = parent1;
            new_X(:,i+1) = parent2;
        end
        cross_sols_remaining = cross_sols_remaining - 2;
    end
    % If odd # of sols needed:
    if cross_sols_remaining == 1

        crossover_decider = rand;
        rand_index1 = randi(n);
        parent1 = mating_pool(:,rand_index1);

        if crossover_decider < cc
            rand_index2 = randi(n);
            parent2 = mating_pool(:,rand_index2);

            for chan = 1:split
                new_X(chan,n) = parent1(chan);
            end
            for chan = split+1:c
                new_X(chan,n) = parent2(chan);
            end
        else
            new_X(:,n) = parent1;
        end
    end

    % Mutation:
    % causes random change of gene based on given probability
    for i = 1:c
        for j = num_elites+1:n
            if rand < mc
                new_X(i,j) = round(rand);
            end
        end
    end

    % Recalculate fitness scores in new_X using fitness function
    for i = 1:n
        fitness_scores(i) = fitness(new_X(:,i),tt,info);
    end
    
    % Update histories of current gen
    history_fscore(:,gen) = fitness_scores;
    history_mpool(:,:,gen) = mating_pool;
    history_fsums(gen) = sum(fitness_scores);

    % Update/Record Matrix X
    X = new_X;
    history_X(:,:,gen) = X;
end