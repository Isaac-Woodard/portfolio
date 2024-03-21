function idx = randi_weight(n, weight)
    % Randomly returns an index n times for one of the elements in weight.
    % Weighting is used to bias the result.
    %
    % Note: The values in weight should add to one.
    arguments (Input)
        n (1,1) int32
        weight (:,1) double
    end
    arguments (Output)
        idx (:,1) int32
    end

    weight_sum = cumsum(weight);
    weight_sum(end) = 1; % Guarantees last value is 1 even if there are floating point errors.
    r = rand([n 1]);

    [~, ~, idx] = histcounts(r, weight_sum);
    idx = idx + 1;
end