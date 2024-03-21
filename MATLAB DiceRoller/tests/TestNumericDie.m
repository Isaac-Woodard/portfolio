classdef TestNumericDie < TestDie
    
    methods(TestClassSetup)
        % Shared setup for the entire test class
    end
    
    methods(TestMethodSetup)
        % Setup for each test
        function create_die(testCase)
            testCase.die = NumericDie([1 2 3 4 5 6]);
        end
    end
    
    methods(Test)
        function test_roll(testCase)
            [result, color] = testCase.die.roll(1);

            testCase.verifyTrue(length(result) == 1);
            testCase.verifyTrue(any(result == testCase.die.Sides))
            testCase.verifyTrue(length(color) == 3);
        end

        function test_roll5(testCase)
            [result, color] = testCase.die.roll(5);

            testCase.verifyTrue(length(result) == 5);
            testCase.verifyTrue(length(result) == length(color))
        end
    end
end