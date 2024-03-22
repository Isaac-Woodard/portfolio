%% Roll a d6 5 times.
d6 = NumericDie(1:6);
results1 = d6.roll(5);
disp(results1)

%% Roll a weather categoric die 5 times.
weather_die = CategoricDie(["clear" "clouds" "rain" "snow"]);
results2 = weather_die.roll(5);
disp(results2)

%% Roll five d6 once. 
d6 = NumericDie(1:6);
dice = repmat(d6, [5 1]);
results3 = zeros(5, 1);

for i = 1:length(dice)
    results3(i) = dice(i).roll();
end

disp(results3)

% Note: Less efficient than rolling the same die multiple times since 
% NumericDie isn't vectorized.

%% Roll a weighted numeric die and plot a histogram of results.
weight_step = [1 1 10 1 1 1 1 1 1 1]';
d10_step = NumericDie(1:10, [0 0 0], weight);
d10_step.roll(100);
histogram(d10_step.History); 
title("Step Weighted Die"); xlabel("Value"); ylabel("Count");

%%  Roll a Gaussian weighted numeric die and plot a histogram of results.
sides = 10;
weight_gaussian = probability.normdist_discrete(sides, 6);
d10_gaussian = NumericDie(1:sides, [0 0 0], weight_gaussian);
d10_gaussian.roll(1000);
histogram(d10_gaussian.History, sides); 
title("Gaussian Weighted Die"); xlabel("Value"); ylabel("Count");