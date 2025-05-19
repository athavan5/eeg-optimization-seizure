tt = edfread('chb01_01.edf');
info = edfinfo('chb01_01.edf');

% Parameters
n = 10; % # of initial solutions
c = 23; % # of channels

num_elites = 2;
num_gens = 1;

% Create random c x n matrix of 0s and 1s (X)
X = round(rand(c,n));
new_X = zeros(size(X));

% Create fitness score vector of size n x 1
fitness_scores = zeros(n,1);

% Repeat algorithm for specified number of generations
for gen = 1:num_gens

    % Calculate fitness scores using fitness function
    for i = 1:n
        fitness_scores(i) = fitness(X(:,i));
    end
    fs = fitness_scores;
    
    % Sort matrix into new matrix based on fitness scores
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
    fitness_scores = fs;
end

