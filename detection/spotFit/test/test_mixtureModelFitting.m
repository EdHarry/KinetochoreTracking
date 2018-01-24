%TEST_MIXTUREMODELFITTING is a test function for mixtureModelFitting using the xUnit framework
%
% To run all tests, call
%    runtests test_mixtureModelFitting
% To run a specific test, call
%    test_mixtureModelFitting(specificTest)
%
% SEE ALSO mixtureModelFitting, xUnit
%
% created with MATLAB ver.: 7.14.0.739 (R2012a) on Mac OS X  Version: 10.6.8 Build: 10K549 
%
% created by: Jacques Boisvert
% DATE: 24-May-2012
%
% Last revision $Rev: 2723 $ $Date: 2012-05-03 03:33:21 -0400 (Thu, 03 May 2012) $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
classdef test_mixtureModelFitting < TestCase
	properties
		sampleProperty %defined by setup method
	end

	methods
%% Constructor
		function obj = test_mixtureModelFitting(name)
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
			mixtureModelFitting
		end
	end %methods
end % classdef% function out = setup
