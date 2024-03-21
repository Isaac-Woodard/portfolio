classdef TestDie < matlab.unittest.TestCase
    
    methods(TestClassSetup)
        % Shared setup for the entire test class.
        
    end
    
    methods(TestMethodSetup)
        % Setup for each test.
    end
    
    methods(Test)
        % Test methods.

        function test_add_to_history(testCase)
            result = testCase.die.roll(1);
            testCase.verifyTrue(result == testCase.die.History);
        end
        
        function test_clear_history(testCase)
            testCase.verifyFail("Unimplemented test");
        end
    end
end