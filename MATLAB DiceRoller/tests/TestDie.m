classdef TestDie < matlab.unittest.TestCase
    
    methods(TestClassSetup)
        % Shared setup for the entire test class.
        
    end
    
    methods(TestMethodSetup)
        % Setup for each test.
    end
    
    methods(Test)
        function test_add_to_history(testCase)
            result = testCase.die.roll(1);

            testCase.verifyEqual(result, testCase.die.History);
        end
        
        function test_clear_history(testCase)
            testCase.die.roll(1);
            testCase.die.clear_history();
            
            testCase.verifyEmpty(testCase.die.History);
        end
    end
end