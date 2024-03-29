classdef TestCategoricDie < TestDie

    methods(TestMethodSetup)
        function create_die(testCase)
            testCase.Die = CategoricDie(["a" "b" "c" "d"]);
        end
    end
    
    methods(Test)
        function test_roll(testCase)
            [result, color] = testCase.Die.roll(1);

            testCase.verifyEqual(length(result), 1);
            testCase.verifyTrue(any(strcmp( ...
                result, string(testCase.Die.Sides))));
            testCase.verifyEqual(length(color), 3);
        end

        function test_roll5(testCase)
            [result, color] = testCase.Die.roll(5);

            testCase.verifyEqual(length(result), 5);
            testCase.verifyEqual(length(result), length(color))
        end
    end
end