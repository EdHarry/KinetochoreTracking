%TEST_CLUSTERSPOT is a test function for clusterSpot using the xUnit framework
%
% To run all tests, call
%    runtests test_clusterSpot
% To run a specific test, call
%    test_clusterSpot(specificTest)
%
% SEE ALSO clusterSpot, xUnit
%
% created with MATLAB ver.: 7.14.0.739 (R2012a) on Mac OS X  Version: 10.6.8 Build: 10K549
%
% created by: Jacques Boisvert
% DATE: 24-May-2012
%
% Last revision $Rev: 2723 $ $Date: 2012-05-03 03:33:21 -0400 (Thu, 03 May 2012) $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
classdef test_clusterSpot < TestCase
    properties
        sampleProperty %defined by setup method
        imgTest1
        imgTest2
        imgTest3
        cl1
        cl2
        cl3
        coord
        coord2
        coord3
    end
    
    methods
        %% Constructor
        function obj = test_clusterSpot(name)
            % Argument 'name' allows running specific tests
            obj = obj@TestCase(name);
        end
        %% Setup
        function setUp(obj)
            fp = which('testClusterFile.mat');
            load(fp);
            obj.imgTest1 = img;
            obj.imgTest2 = img2;
            obj.imgTest3 = img3;
            obj.cl1 = cl; obj.cl2 = cl2; obj.cl3 = cl3;
            obj.coord = coord; obj.coord2 = coord2 ; obj.coord3 = coord3;
            % common setup code goes here
        end
        %% Teardown
        function tearDown(obj)
            % cleanup code goes here
        end
        %% Test #1
        %Test with an image consisting of 2 clusters of 18 spots each.
        function testNumberOne(obj)
            % first test. All test methods need to start with 'test'
            %Load image test -> clusterTestImage.mat ->Contain the test
            %image under the name img.
            %Spots were made with a gaussian sigma of 2.
            %Test file also has the coordinate and the expected answer.
            mp = distinguishable_colors(size(unique(obj.cl1),1));
            dimshow(obj.imgTest1);hold on; scatter(obj.coord(:,2),obj.coord(:,1),[],mp(obj.cl1,:));
            c = clusterSpot(obj.coord,2);
            mp = distinguishable_colors(size(unique(c),1));
            dimshow(obj.imgTest1);hold on; scatter(obj.coord(:,2),obj.coord(:,1),[],mp(c,:));
            assertEqual(c,obj.cl1);
        end
        
        %% Test # 2
        function testNumberTwo(obj)
            mp = distinguishable_colors(size(unique(obj.cl2),1));
            dimshow(obj.imgTest2);hold on; scatter(obj.coord2(:,2),obj.coord2(:,1),[],mp(obj.cl2,:));
            c = clusterSpot(obj.coord2,2);
            mp = distinguishable_colors(size(unique(c),1));
            dimshow(obj.imgTest2);hold on; scatter(obj.coord2(:,2),obj.coord2(:,1),[],mp(c,:));
            assertEqual(c,obj.cl2);
        end
        
        function testNumberThree(obj)
            mp = distinguishable_colors(size(unique(obj.cl3),1));
            dimshow(obj.imgTest3);hold on; scatter(obj.coord3(:,2),obj.coord3(:,1),[],mp(obj.cl3,:));
            c = clusterSpot(obj.coord3,2);
            mp = distinguishable_colors(size(unique(c),1));
            dimshow(obj.imgTest3);hold on; scatter(obj.coord3(:,2),obj.coord3(:,1),[],mp(c,:));
            assertEqual(c,obj.cl3);
        end
    end %methods
end % classdef% function out = setup
