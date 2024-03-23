classdef TestDie < matlab.unittest.TestCase

    properties
        Die
    end

    methods(TestClassSetup)
        function updir(testCase)
            cd("..\");
        end
    end

    methods(TestClassTeardown)
        function downdir(testCase)
            cd("tests");
        end
    end
    
    methods(TestMethodSetup, Abstract)
        create_die(testCase)
    end
    
    methods(Test)
        function test_add_to_history(testCase)
            result = testCase.Die.roll(1);

            testCase.verifyEqual(result, testCase.Die.History);
        end
        
        function test_clear_history(testCase)
            testCase.Die.roll(1);
            testCase.Die.clear_history();
            
            testCase.verifyEmpty(testCase.Die.History);
        end
    end
end