classdef TestCategoricDie < TestDie
    
    methods(TestClassSetup)
        % Shared setup for the entire test class
    end
    
    methods(TestMethodSetup)
        % Setup for each test
    end
    
    methods(Test)
        % Test methods
        
        function test_roll(testCase)
            testCase.verifyFail("Unimplemented test");
        end

        function test_add_to_history(testCase)
            %NOTE: The test wil have to be indirect since add_to_history() is private.
            testCase.verifyFail("Unimplemented test");
        end
    end
end