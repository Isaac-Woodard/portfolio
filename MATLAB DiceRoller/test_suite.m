import matlab.unittest.TestSuite

switch ("All")
    case "All"
        tests = TestSuite.fromFolder("tests");
    case "Categoric"
        tests = TestSuite.fromFile("tests\TestCategoricDie.m");
    case "Numeric"
        tests = TestSuite.fromFile("tests\TestNumericDie.m");
    case "probability"
        tests = TestSuite.fromFile("tests\TestProbability.m");
end

run(tests);