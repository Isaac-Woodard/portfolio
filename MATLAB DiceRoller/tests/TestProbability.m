% Tests functions in the probability package.

classdef TestProbability < matlab.unittest.TestCase
    
    methods(Test)
        function test_normdist_discrete(testCase)
            n = 10; range = 1;
            weights = probability.normdist_discrete(n, range);

            testCase.verifyEqual(length(weights), n);
            testCase.verifyEqual(weights(1), weights(end));
        end

        function test_randi_weight(testCase)
            weight = [0.25 0.25 0.25 0.25];
            idx = probability.randi_weight(10, weight);

            testCase.verifyEqual(length(idx), 10);
            testCase.verifyTrue(max(idx) <= length(weight));
            testCase.verifyTrue(min(idx) > 0)
        end
    end
    
end