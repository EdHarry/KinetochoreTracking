%TEST_TESTDISTANCESIGNIFICANCE is a test function for testDistanceSignificance using the xUnit framework
%
% To run all tests, call
%    runtests test_testDistanceSignificance
% To run a specific test, call
%    test_testDistanceSignificance(specificTest)
%
% SEE ALSO testDistanceSignificance, xUnit
%
% created with MATLAB ver.: 7.14.0.739 (R2012a) on Mac OS X  Version: 10.6.8 Build: 10K549 
%
% created by: Jacques Boisvert
% DATE: 07-Jun-2012
%
% Last revision $Rev: 2842 $ $Date: 2012-06-07 14:20:08 -0400 (Thu, 07 Jun 2012) $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
classdef test_testDistanceSignificance < TestCase
	properties
		sampleProperty %defined by setup method
	end

	methods
%% Constructor
		function obj = test_testDistanceSignificance(name)
		% Argument 'name' allows running specific tests
			obj = obj@TestCase(name);
		end
%% Setup
		function setUp(obj)
			% common setup code goes here
		end
%% Teardown
		function tearDown(obj)
			% cleanup code goes here
		end
%% Test #1
		function testNumberOne(obj)
			% first test. All test methods need to start with 'test'
			testDistanceSignificance
		end
	end %methods
end % classdef% function out = setup
