function weights = normdist_discrete(n, range, options)
    % Returns a vector of discrete weights based on the noraml distribution.
    % The weights are centered about mean=0 with standard deviation=1.
    arguments (Input)
        n (1,1) int32
        range (1,1) double {mustBePositive}
        options.mean (1,1) double = 0
        options.sigma (1,1) double {mustBePositive} = 1
    end
    arguments (Output)
        weights (:,1) double
    end

    x = linspace(-0.5*range, 0.5*range, n);

    a = (1 / (options.sigma * (2 * pi).^0.5));
    b = -0.5 * ((x - options.mean) / options.sigma).^2;
    weights = a * exp(b);
end