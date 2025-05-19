function fitness_score = fitness(x)
    x_size = size(x);
    total_channels = x_size(1);
    
    % Calculate Score based on # of Channels Used (from 0-1; less is better)
    channels_used = sum(x);
    chan_score = 1 - channels_used/total_channels;

    % Calculate Score based on Accuracy (tbd)
    acc_score = 0;

    fitness_score = chan_score + acc_score;
end