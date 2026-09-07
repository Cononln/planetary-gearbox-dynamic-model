function tests = test_stage1
tests = functiontests(localfunctions);
end

function testValidationRuns(testCase)
report = step1_validate_18dof;
verifyEqual(testCase,report.nDof,18);
verifyLessThan(testCase,report.periodicResidual,1e-12);
end
