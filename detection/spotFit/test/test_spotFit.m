%TEST_SPOTFIT is a test function for spotFit using the xUnit framework
%
% To run all tests, call
%    runtests test_spotFit
% To run a specific test, call
%    test_spotFit(specificTest)
%
% SEE ALSO spotFit, xUnit
%
% created with MATLAB ver.: 7.14.0.739 (R2012a) on Mac OS X  Version: 10.6.8 Build: 10K549
%
% created by: Jacques Boisvert
% DATE: 24-May-2012
%
% Last revision $Rev: 2723 $ $Date: 2012-05-03 03:33:21 -0400 (Thu, 03 May 2012) $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
classdef test_spotFit < TestCase
    properties
        sampleProperty %defined by setup method
    end
    
    methods
        %% Constructor
        function obj = test_spotFit(name)
            % Argument 'name' allows running specific tests
            obj = obj@TestCase(name);
        end
        %% Setup
        function setUp(obj)
            % common setup code goes here
            %Set seed to have reproductible test case
            seed = RandStream('mt19937ar','seed',1);
            RandStream.setGlobalStream(seed);
        end
        %% Teardown
        function tearDown(obj)
            % cleanup code goes here
        end
        %% Test One Spot One Kernel Perfect Start
        function testOneSpotOneKernel(obj)
            % first test. All test methods need to start with 'test'
            
            %nTest test on pixel locking.
            nTest = 10^3;
            %2D
            coordList = zeros(nTest,2);
            groundTruthCoordList = zeros(nTest,2);
            fullMatrixCoordMean = zeros(1,2);
            %
            %3D
            %             coordList = zeros(nTest,3);
            %             groundTruthCoordList = zeros(nTest,3);
            %             fullMatrixCoordMean = zeros(1,3);
            for cTest = 1:nTest
                %create test image
                fprintf('Run #%d\n',cTest);
                %2D
                imgSize = [64 64];
                coordinate = [ imgSize(1)/2,imgSize(2)/2];
                coordinate = coordinate + (rand(1,2)-0.5)*6;
                
                %3D
                %                 %Image size
                %                 imgSize = [ 64 64 64];
                %                 %Starting coordinate is in the middle
                %                 coordinate = [ imgSize(1)/2, imgSize(2)/2, imgSize(3)/2];
                %                 %Add variation in the center position
                %                coordinate = coordinate + (rand(1,3)-0.5)*6;
                
                
                %Save ground truth
                groundTruthCoordList(cTest,:) = coordinate;
                %Signal to noise Ratio
                snr = 10^3;
                %Spot size
                featureSize = 2;
                
                [nseImg,~,~] = createMixtureFitTestImage('coordinates',coordinate,...
                    'snr',snr,'imgSize',imgSize,'featureSize',featureSize);
                %Maddox
                [newCoord,fullAmp,fullBg,fullMatrix] = spotMMFit(nseImg,coordinate,featureSize,'fitNPlusOne',1);
                %LCCB
                %coordinate = coordinate - size(nseImg)/2;
                %                coordinate = cat(2,coordinate,[1,2,0]);
                %                [prm,prmStd,C] = fitGaussianMixture2D(nseImg,coordinate,'xyac');
                %                fullMatrix = C;
                %                newCoord = [prm(:,1) prm(:,2)];
                %                fullAmp = prm(:,3);
                %                fullBg = prm(:,4);
                % LCCB END
                
                fullMatrixCoordMean(1) = fullMatrixCoordMean(1) + fullMatrix(1,1);
                fullMatrixCoordMean(2) = fullMatrixCoordMean(2) + fullMatrix(2,2);
                %3D
                %                fullMatrixCoordMean(3) = fullMatrixCoordMean(3) + fullMatrix(3,3);
                
                %Save output coordinates
                if ~(isempty(newCoord))
                    coordList(cTest,:) = newCoord;
                else
                    coordList(cTest,:) = NaN;
                end
                
            end
            
            fullMatrixCoordMean = fullMatrixCoordMean ./ nTest;
            
            %2D
            %Plot the residual of the coordinates.
            dfig,plot( ( coordList - groundTruthCoordList) ) ;
            %Normalized
            coordList(:) = coordList(:) - floor(coordList(:));
            groundTruthCoordList(:) = groundTruthCoordList(:) - floor(groundTruthCoordList(:));
            %Plot found coordinates
            dfig,plot(coordList(:,1),coordList(:,2),'ob');
            hold on
            %             Overlap with ground truth coordinate
            plot(groundTruthCoordList(:,1),groundTruthCoordList(:,2),'xg');
            axis equal;
            %3D
            %             dfig;plot((coordList-groundTruthCoordList));
            %             coordList(:) = coordList(:) - floor(coordList(:));
            %             groundTruthCoordList(:) = groundTruthCoordList(:) - floor(groundTruthCoordList(:));
            %             dfig;plot3(coordList(:,1),coordList(:,2),coordList(:,3),'ob');
            %             hold on; plot3(groundTruthCoordList(:,1),groundTruthCoordList(:,2),groundTruthCoordList(:,3),'xg');
            %             axis equal;
            
        end
        
        %% test One Spot One Kernel Off Start
        function testOneSpotOneKernelOffStart(obj)
            %n test on offStart and Signal-to-Noise ratio.
            nTest = 10^2;
            
            dataIdx = 1;
            %2D
            data = zeros(size((1.6:.2:10),2),9);
            %3D
            %data = zeros(size((1.6:.2:10),2),11);
            for snr = 1.6:.2:10
                fprintf('SNR : %.3f\n',snr);
                
                %coordList ->output coord
                %2D
                fullMatrixMean = zeros(1,4);
                coordList = NaN(nTest,2);
                inputCoordList = NaN(nTest,2);
                groundTruthCoordList = NaN(nTest,2);
                goneCoord = NaN(nTest,2);
                
                %3D
                %                 coordList = NaN(nTest,3);
                %                 fullMatrixMean = zeros(1,5);
                %                 inputCoordList = NaN(nTest,3);
                %                 groundTruthCoordList = NaN(nTest,3);
                %                 goneCoord = NaN(nTest,3);
                
                ampList = NaN(nTest,1);
                bgList = NaN(nTest,1);
                for cTest = 1:nTest
                    
                    %create test image
                    %2D
                    imgSize = [64 64];
                    coordinate = [32 32];
                    
                    %3D
                    % imgSize = [64 64 64];
                    %coordinate = [32 32 32];
                    
                    %2D
                    offSet = (rand(1,2) - [0.5,0.5])*8;
                    %3D
                    %offSet = (rand(1,3) - [0.5,0.5,0.5])*6;
                    
                    groundTruthCoordList(cTest,:) = coordinate;
                    featureSize = 2;
                    [nseImg,~,~] = createMixtureFitTestImage('coordinates',coordinate,...
                        'snr',snr,'imgSize',imgSize,'featureSize',featureSize);
                    
                    
                    inputCoordList(cTest,:) = (coordinate+offSet);
                    %MADDOX
                    [newCoord,amp,bg,fullMatrix] = spotMMFit(nseImg,(coordinate+offSet),featureSize,'fitNPlusOne',1,'',snr);
                    
                    %LCCB
                    %                     inputCoordinate = cat(2,coordinate+offSet,[1,2,0]);
                    %                     inputCoordinate(1,1:2) = inputCoordinate(1,1:2) - size(nseImg)/2;
                    %                     [prm,~,C] = fitGaussianMixture2D(nseImg,inputCoordinate,'xyac');
                    %                     fullMatrix = C;
                    %                     newCoord = [prm(:,1),prm(:,2)];
                    %                     amp = prm(:,3);
                    %                     bg = prm(:,4);
                    %                     newCoord = newCoord + size(nseImg)/2;
                    %LCCB END
                    
                    if ~(isempty(newCoord))
                        %2D
                        fullMatrixMean(1) = fullMatrixMean(1) + fullMatrix(1,1);
                        fullMatrixMean(2) = fullMatrixMean(2) + fullMatrix(2,2);
                        fullMatrixMean(3) = fullMatrixMean(3) + fullMatrix(3,3);
                        fullMatrixMean(4) = fullMatrixMean(4) + fullMatrix(4,4);
                        
                        %3D
                        %fullMatrixMean(5) = fullMatrixMean(5) + fullMatrix(5,5);
                        
                        if size(newCoord,1) > 1
                            disp('Nooo!');
                        end
                        
                        coordList(cTest,:) = newCoord;
                        ampList(cTest,1) = amp;
                        bgList(cTest,1) = bg;
                        
                    else
                        groundTruthCoordList(cTest,:) = NaN;
                        goneCoord(cTest,:) = inputCoordList(cTest,:);
                        coordList(cTest,:) = NaN;
                        ampList(cTest,1) = NaN;
                    end
                end
                
                fullMatrixMean = fullMatrixMean ./nTest;
                
                %2D
                dfig,plot(coordList(:,1),coordList(:,2),'ob');
                title((snr));
                hold on
                plot(groundTruthCoordList(:,1),groundTruthCoordList(:,2),'xg');
                plot(inputCoordList(:,1),inputCoordList(:,2),'or');
                
                plot(goneCoord(:,1),goneCoord(:,2),'xb');
                % axis equal;
                %dfig,plot( ( coordList - groundTruthCoordList) ) ;
                %3D
                %                 dfig,plot3(coordList(:,1),coordList(:,2),coordList(:,3),'ob');
                %                 title(snr);
                %                 hold on
                %                 plot3(groundTruthCoordList(:,1),groundTruthCoordList(:,2),groundTruthCoordList(:,3),'xg');
                %                 plot3(inputCoordList(:,1),inputCoordList(:,2),inputCoordList(:,3),'or');
                %                 plot3(goneCoord(:,1),goneCoord(:,2),goneCoord(:,3),'xb');
                %                 axis equal;
                
                
                
                std = sqrt( nanmean(( (coordList-groundTruthCoordList).^2)) ...
                    - nanmean( (coordList-groundTruthCoordList)).^2    );
                ampStd = sqrt(nanmean((ampList - 1).^2) - nanmean(ampList - 1).^2);
                bgStd = sqrt(nanmean((bgList - 0).^2) - nanmean(bgList -0).^2);
                %2D
                data(dataIdx,1) = snr;
                data(dataIdx,2:3) = std;
                data(dataIdx,4) = ampStd;
                data(dataIdx,6:9) = fullMatrixMean;
                data(dataIdx,5) = bgStd;
                dataIdx = dataIdx + 1;
                
                %3D
                %                 data(dataIdx,1) = snr;
                %                 data(dataIdx,2:4) = std;
                %                 data(dataIdx,5) = ampStd;
                %                 data(dataIdx,6) = bgStd;
                %                 data(dataIdx,7:11) = fullMatrixMean;
                %                dataIdx = dataIdx + 1;
            end
            %2D
            dfig;plot(data(:,1),data(:,2),'ob');
            hold on
            plot(data(:,1),sqrt(data(:,6)),'or');
            hold off
            
            dfig;plot(data(:,1),data(:,3),'ob');
            hold on
            plot(data(:,1),sqrt(data(:,7)),'or');
            
            
            dfig;plot(data(:,1),data(:,4),'ob');
            hold on
            plot(data(:,1),sqrt(data(:,8)),'or');
            
            %3D
            %              dfig; plot(data(:,1),data(:,2),'ob');
            %              hold on
            %              plot(data(:,1),sqrt(data(:,7)),'or');
            %              title('std in X');
            %              dfig;plot(data(:,1),data(:,3),'ob');
            %              hold on
            %              plot(data(:,1),sqrt(data(:,8)),'or');
            %              title('std in Y');
            %              dfig;plot(data(:,1),data(:,4),'ob');
            %              hold on
            %              plot(data(:,1),sqrt(data(:,9)),'or');
            %              title('std in Z');
            %              dfig;plot(data(:,1),data(:,5),'ob');
            %              hold on
            %              plot(data(:,1),sqrt(data(:,10)),'or');
            %              dfig;plot(data(:,1),data(:,6),'ob');
            %              hold on
            %              plot(data(:,1),sqrt(data(:,11)),'or');
            
        end
        
        
        %% test One Spot Two Kernel
        function testOneSpotTwoKernel(obj)
            
            nTest = 50;
            
            dataIdx = 1;
            nDistance = size((0.5:0.5:10),2);
            %nDistance = size((4:0.5:10),2);
            nSnr = size((2:.2:10),2);
            %nSnr = size((8:.2:10),2);
            data = zeros(nSnr+nDistance,7);
            
            for snr = 2:.2:10
                fprintf('Snr test : %.2f\n',snr);
                
                
                for distance = 0.5:0.5:10
                    fprintf('Distance test : %.2f\n',distance);
                    nHit1 = 0;
                    nHit2 = 0;
                    nHit0 = 0;
                    
                    groundTruthCoordList = zeros(nTest,2);
                    coordList = NaN(nTest,2);
                    for cTest = 1:nTest
                        %create test image
                        imgSize = [ 64 64];
                        coordinate = [ 32 32];
                        groundTruthCoordList(cTest,:) = coordinate;
                        featureSize = 2;
                        [nseImg,~,~] = createMixtureFitTestImage('coordinates',coordinate,...
                            'snr',snr,'imgSize',imgSize,'featureSize',featureSize);
                        inputCoord = [32,32;32+distance,32];
                        %MADDOX
                        [newCoord,newAmp,newBg,fullMatrix] = spotMMFit(nseImg,inputCoord,featureSize,...
                            'fitNPlusOne',1,'amplitudeCutOff',0.01);
                        
                        %LCCB
                        %                         inputCoordinate = cat(2,inputCoordinate,[1,0,0;1,0,0]);
                        %                         [prm,prmStd,C] = fitGaussianMixture2D(nseImg,inputCoordinate,'xy');
                        %                         C = fullMatrix;
                        %                         newAmp = prm(:,3);
                        %                         newBg = prm(:,4);
                        %                         newCoord = [prm(:,1),prm(:,2)];
                        %LCCB END
                        
                        if size(newCoord,1) == 2
                            nHit2 = nHit2+1;
                        elseif size(newCoord,1) == 1
                            nHit1 = nHit1+1;
                            coordList(cTest,:) = newCoord;
                        elseif size(newCoord,1) == 0
                            nHit0 = nHit0+1;
                        end
                    end
                    
                    nHit0 = nHit0/nTest;
                    nHit1 = nHit1/nTest;
                    nHit2 = nHit2/nTest;
                    
                    std = sqrt( nanmean(( (coordList-groundTruthCoordList).^2)) - nanmean( (coordList-groundTruthCoordList)).^2    );
                    
                    data(dataIdx,1) = snr;
                    data(dataIdx,2) = distance;
                    data(dataIdx,3:4) = std;
                    data(dataIdx,5) = nHit0;
                    data(dataIdx,6) = nHit1;
                    data(dataIdx,7) = nHit2;
                    dataIdx = dataIdx + 1;
                    
                    %                      dfig;plot(coordList(:,1),coordList(:,2),'og');
                    %                      hold on;
                    %                     plot(groundTruthCoordList(:,1),groundTruthCoordList(:,2),'xr');
                end
            end
            z = reshape(data(:,6),nDistance,nSnr);
            x = 2:.2:10;
            y = 0.5:0.5:10;
            dfig,contourf(x,y,z);
            title('% of test with 1 kernel left');
        end
        
        %% test One Spot One Kernel Wrong Sigma
        function testOneSpotOneKernelWrongSigma(obj)
            nTest = 100;
            data = zeros(size(1.6:.2:10,2)+ size(1:.2:6,2),6);
            dataIdx = 1;
            nSnr = size((1.6:.2:10),2);
            nSigma = size((1:.2:6),2);
            for snr = 1.6:.2:10
                fprintf('snr test : %.2f\n',snr);
                for sigma = 1:.2:6
                    fprintf('sigma : %.2f\n',sigma);
                    
                    coordList = NaN(nTest,2);
                    groundTruthCoordList = zeros(nTest,2);
                    for cTest = 1 :nTest
                        %create test image
                        imgSize = [ 64 64];
                        coordinate = [ 32 32];
                        groundTruthCoordList(cTest,:) = coordinate;
                        featureSize = 2;
                        [nseImg,~,~] = createMixtureFitTestImage('coordinates',coordinate,...
                            'snr',snr,'imgSize',imgSize,'featureSize',featureSize);
                        [newCoord,amp,bg,fullMatrix] = spotMMFit(nseImg,coordinate,sigma,'fitNPlusOne',1);
                        if ~isempty(newCoord)
                            if size(newCoord,1) > 1
                                fprintf('%d %d',newCoord(1),newCoord(2));
                            else
                                coordList(cTest,:) = newCoord;
                            end
                        end
                    end
                    std = sqrt( nanmean(( (coordList-groundTruthCoordList).^2)) - nanmean( (coordList-groundTruthCoordList)).^2 );
                    data(dataIdx,1) = snr;
                    data(dataIdx,2) = sigma;
                    data(dataIdx,3:4) = std;
                    dataIdx = dataIdx + 1;
                    
                    %                     dfig;scatter(groundTruthCoordList(:,2),groundTruthCoordList(:,1),'x');
                    %                     hold on
                    %                     scatter(coordList(:,2),coordList(:,1),'o');
                    %                     title(sprintf('snr : %.2f , sigma : %.2f',snr,sigma));
                end
            end
            xStd = reshape(data(:,3),nSigma,nSnr);
            yStd = reshape(data(:,4),nSigma,nSnr);
            for col = 1:nSnr
                dfig;
                subplot(1,2,1),plot(xStd(:,col),'o');
                subplot(1,2,2), plot(yStd(:,col),'o');
            end
        end
        
        
        
        
        %% test Two Spots Two Kernels Perfect Start
        function testTwoSpotTwoKernel(obj)
            nTest = 100;
            %nSnr = size((1.6:.2:6),2);
            nSnr = size((8:.2:10),2);
            nDistance = size((4:0.5:10),2);
            %nDistance = size((1:.1:10),2);
            data = zeros(nSnr + nDistance,9);
            dataIdx = 1;
            
            for snr = 8:.2:10
                fprintf('snr test : %.2f\n',snr);
                
                for distance = 4:0.5:10
                    
                    fprintf('distance : %.2f\n',distance);
                    groundTruthCoordList = zeros(nTest,4);
                    coordList = zeros(nTest,4);
                    nHit0 = 0;
                    nHit1 = 0;
                    nHit2 = 0;
                    for cTest = 1:nTest
                        imgSize = [ 64 64];
                        coordinate = [ 22 22;22+distance 22];
                        groundTruthCoordList(cTest,1:2) = coordinate(1,:);
                        groundTruthCoordList(cTest,3:4) = coordinate(2,:);
                        featureSize = 2;
                        [nseImg,~,~] = createMixtureFitTestImage('coordinates',coordinate,...
                            'snr',snr,'imgSize',imgSize,'featureSize',featureSize);
                        %MADDOX
                        [newCoord,amp,bg,fullMatrix] = spotMMFit(nseImg,coordinate,2,'fitNPlusOne',1);
                        %LCCB
                        %                         inputCoordinate = cat(2,coordinate(1,:),1);
                        %                         inputCoordinate = cat(2,inputCoordinate,coordinate(2,:),[1 2 0]);
                        %                         [prm,~,C] = fitGaussianMixture2D(nseImg,inputCoordinate,'xyac');
                        %                         fullMatrix = C;
                        %                         newCoord = [prm(:,1),prm(:,2)];
                        %                         amp = prm(:,3);
                        %                         bg = prm(:,4);
                        %LCCB END
                        if size(newCoord,1) == 2
                            coordList(cTest,1:2) = newCoord(1,:);
                            coordList(cTest,3:4) = newCoord(2,:);
                            nHit2 = nHit2+1;
                        end
                        if size(newCoord,1) == 1
                            nHit1 = nHit1+1;
                        end
                        if size(newCoord,1) == 0
                            nHit0 = nHit0+1;
                        end
                    end
                    nHit0 = nHit0/nTest;
                    nHit1 = nHit1/nTest;
                    nHit2 = nHit2/nTest;
                    
                    std = sqrt( nanmean(( (coordList-groundTruthCoordList).^2)) - nanmean( (coordList-groundTruthCoordList)).^2    );
                    dfig;
                    plot(coordList(:,1),coordList(:,2),'ob');
                    hold on
                    plot(coordList(:,3),coordList(:,4),'or');
                    plot(groundTruthCoordList(:,1),groundTruthCoordList(:,2),'xg');
                    plot(groundTruthCoordList(:,3),groundTruthCoordList(:,4),'xg');
                    axis equal;
                    data(dataIdx,3:6) = std;
                    data(dataIdx,1) = snr;
                    data(dataIdx,7) = nHit0;
                    data(dataIdx,8) = nHit1;
                    data(dataIdx,9) = nHit2;
                    data(dataIdx,2) = distance;
                    dataIdx = dataIdx +1;
                    
                    
                end
            end
            dfig;
            hit2 = reshape(data(:,end),nDistance,nSnr);
            hit1 = reshape(data(:,end-1),nDistance,nSnr);
            contourf(hit2);
            title('nHit = 2');
            
            dfig; contourf(hit1);
            title('nHit = 1');
            
        end
        
        %% test Two Kernel One Spot with different sigma
        function testTwoSpotOneKernelSigma(obj)
            
            nTest = 100;
            
            nSnr = size((2:.2:10),2);
            nSigma = size((2:.2:10),2);
            data = zeros(nSnr+nSigma,9);
            dataIdx = 1;
            for snr = 2:.2:10
                for featureSize = 2:.2:10
                    fprintf('snr test : %.2f\n',snr);
                    
                    coordList = NaN(nTest,4);
                    groundTruthCoordList = zeros(nTest,4);
                    centroidList = NaN(nTest,4);
                    nHit2 = 0;
                    nHit1 = 0;
                    nHit3 = 0;
                    for cTest = 1 :nTest
                        %create test image
                        imgSize = [ 128 128];
                        coordinate = [ 64 64;64 64];
                        coordinate(1,:) = coordinate(1,:) - 1.5*featureSize;
                        coordinate(2,:) = coordinate(2,:) + 1.5*featureSize;
                        groundTruthCoordList(cTest,1:2) = coordinate(1,:);
                        groundTruthCoordList(cTest,3:4) = coordinate(2,:);

                        [nseImg,~,~] = createMixtureFitTestImage('coordinates',coordinate,...
                            'snr',snr,'imgSize',imgSize,'featureSize',featureSize);
                        [newCoord,~,~,fullMatrix,newCentroidList] = spotMMFit(nseImg,coordinate(1,:),featureSize,'fitNPlusOne',1,'debug',1);
                        if ~isempty(newCoord)
                            if size(newCoord,1) == 2
                                coordList(cTest,1:2) = newCoord(1,:);
                                coordList(cTest,3:4) = newCoord(2,:);
                                nHit2 = nHit2 + 1;
                                centroidList(cTest,3:4) = newCentroidList(1,:);
                            elseif size(newCoord,1) == 1
                                coordList(cTest,1:2) = newCoord;
                                nHit1 = nHit1 + 1;
                                centroidList(cTest,1:2) = newCentroidList(1,:);
                            end
                            if(size(newCoord,1) > 2)
                                nHit3 = nHit3 + 1;
                                fprintf('more then 2 spots');
                            end
                        end
                    end
                    nHit2 = nHit2 / nTest;
                    nHit1 = nHit1 / nTest;
                    nHit3 = nHit3 / nTest;
                    %dfig; plot(groundTruthCoordList(:,1),groundTruthCoordList(:,2),'xg');
                    %hold on
                    %plot(groundTruthCoordList(:,3),groundTruthCoordList(:,4),'xg');
                    %plot(coordList(:,1),coordList(:,2),'ob')
                    %plot(coordList(:,3),coordList(:,4),'or');
                    %title(snr);
                    centroidList = centroidList + randn(size(centroidList)) .* 0.05;
                    %plot(centroidList(:,1),centroidList(:,2),'.r');
                    %plot(centroidList(:,3),centroidList(:,4),'.g');
                    data(dataIdx,1) = snr;
                    data(dataIdx,2) = nHit1;
                    data(dataIdx,3) = nHit2;
                    data(dataIdx,8) = nHit3;
                    data(dataIdx,9) = featureSize;
                    dataIdx = dataIdx + 1;
                end
            end
            dfig;plot(data(:,1),data(:,2));
            title('% of 1 hit');
            dfig,plot(data(:,1),data(:,3));
            title('% of 2 hit');
            dfig,plot(data(:,1),data(:,8));
            title('% of 3+ hit');
        end
        
        function testTwoSpotOneKernel(obj)
            
            nTest = 1000;
            
            nSnr = size((2:.2:10),2);
            data = zeros(nSnr,8);
            dataIdx = 1;
            for snr = 2:.2:10
                fprintf('snr test : %.2f\n',snr);
                
                coordList = NaN(nTest,4);
                groundTruthCoordList = zeros(nTest,4);
                centroidList = NaN(nTest,4);
                nHit2 = 0;
                nHit1 = 0;
                nHit3 = 0;
                for cTest = 1 :nTest
                    %create test image
                    imgSize = [ 64 64];
                    coordinate = [ 32 32;38 38];
                    groundTruthCoordList(cTest,1:2) = coordinate(1,:);
                    groundTruthCoordList(cTest,3:4) = coordinate(2,:);
                    featureSize = 2;
                    [nseImg,~,~] = createMixtureFitTestImage('coordinates',coordinate,...
                        'snr',snr,'imgSize',imgSize,'featureSize',featureSize);
                    [newCoord,~,~,fullMatrix,newCentroidList] = spotMMFit(nseImg,coordinate(1,:),2,'fitNPlusOne',1,'debug',1);
                    if ~isempty(newCoord)
                        if size(newCoord,1) == 2
                            coordList(cTest,1:2) = newCoord(1,:);
                            coordList(cTest,3:4) = newCoord(2,:);
                            nHit2 = nHit2 + 1;
                            centroidList(cTest,3:4) = newCentroidList(1,:);
                        elseif size(newCoord,1) == 1
                            coordList(cTest,1:2) = newCoord;
                            nHit1 = nHit1 + 1;
                            centroidList(cTest,1:2) = newCentroidList(1,:);
                        end
                        if(size(newCoord,1) > 2)
                            nHit3 = nHit3 + 1;
                            fprintf('more then 2 spots');
                        end
                    end
                end
                nHit2 = nHit2 / nTest;
                nHit1 = nHit1 / nTest;
                nHit3 = nHit3 / nTest;
                dfig; plot(groundTruthCoordList(:,1),groundTruthCoordList(:,2),'xg');
                hold on
                plot(groundTruthCoordList(:,3),groundTruthCoordList(:,4),'xg');
                plot(coordList(:,1),coordList(:,2),'ob')
                plot(coordList(:,3),coordList(:,4),'or');
                title(snr);
                centroidList = centroidList + randn(size(centroidList)) .* 0.05;
                plot(centroidList(:,1),centroidList(:,2),'.r');
                plot(centroidList(:,3),centroidList(:,4),'.g');
                data(dataIdx,1) = snr;
                data(dataIdx,2) = nHit1;
                data(dataIdx,3) = nHit2;
                data(dataIdx,8) = nHit3;
                dataIdx = dataIdx + 1;
            end
            dfig;plot(data(:,1),data(:,2));
            title('% of 1 hit');
            dfig,plot(data(:,1),data(:,3));
            title('% of 2 hit');
            dfig,plot(data(:,1),data(:,8));
            title('% of 3+ hit');
        end
        
        
        
        
        %% test Two Spots Three Kernel
        function testTwoSpotThreeKernel(obj)
            
            nTest = 20;
            
            dataIdx = 1;
            nDistance = size((1:0.5:10),2);
            nSnr = size((2:.2:10),2);
            
            
            data = zeros(nDistance + nSnr,10);
            for snr = 2:.2:10
                fprintf('Snr test : %.2f\n',snr);
                
                
                for distance = 1:0.5:10
                    fprintf('Distance test : %.2f\n',distance);
                    nHit1 = 0;
                    nHit2 = 0;
                    nHit0 = 0;
                    nHit3 = 0;
                    groundTruthCoordList = zeros(nTest,4);
                    coordList = NaN(nTest,6);
                    for cTest = 1:nTest
                        %create test image
                        imgSize = [ 64 64];
                        coordinate = [ 24 24; 32 32];
                        groundTruthCoordList(cTest,1:2) = coordinate(1,:);
                        groundTruthCoordList(cTest,3:4) = coordinate(2,:);
                        featureSize = 2;
                        [nseImg,~,~] = createMixtureFitTestImage('coordinates',coordinate,...
                            'snr',snr,'imgSize',imgSize,'featureSize',featureSize);
                        inputCoord = [26,26;32,32;32+distance,32];
                        
                        [newCoord,~,~,fullMatrix] = spotMMFit(nseImg,inputCoord,featureSize,'fitNPlusOne',1);
                        if ~isempty(newCoord)
                            if size(newCoord,1) == 2
                                nHit2 = nHit2+1;
                                coordList(cTest,1:2) = newCoord(1,:);
                                coordList(cTest,3:4) = newCoord(2,:);
                            elseif size(newCoord,1) == 1
                                coordList(cTest,1:2) = newCoord(1,:);
                                nHit1 = nHit1+1;
                            elseif size(newCoord,1) == 0
                                nHit0 = nHit0+1;
                            elseif size(newCoord,1) == 3
                                nHit3 = nHit3+1;
                                coordList(cTest,1:2) = newCoord(1,:);
                                coordList(cTest,3:4) = newCoord(2,:);
                                coordList(cTest,5:6) = newCoord(3,:);
                            end
                        end
                    end
                    nHit3 = nHit3/nTest;
                    nHit0 = nHit0/nTest;
                    nHit1 = nHit1/nTest;
                    nHit2 = nHit2/nTest;
                    
                    std = sqrt( nanmean(( (coordList(:,1:4)-groundTruthCoordList).^2)) ...
                        - nanmean( (coordList(:,1:4)-groundTruthCoordList)).^2);
                    
                    
                    data(dataIdx,1) = snr;
                    data(dataIdx,2) = distance;
                    data(dataIdx,3:6) = std;
                    data(dataIdx,7) = nHit0;
                    data(dataIdx,8) = nHit1;
                    data(dataIdx,9) = nHit2;
                    data(dataIdx,10) = nHit3;
                    dataIdx = dataIdx + 1;
                    
                    dfig;plot(coordList(:,1),coordList(:,2),'ob');
                    hold on
                    plot(groundTruthCoordList(:,1),groundTruthCoordList(:,2),'xg');
                    plot(groundTruthCoordList(:,3),groundTruthCoordList(:,4),'xg');
                    plot(coordList(:,3),coordList(:,4),'or');
                    plot(coordList(:,5),coordList(:,6),'ok');
                    axis equal;
                end
                %dfig;plot(data(data(:,1) == snr,2),data(data(:,1) == snr,9));
                %title(snr);
            end
            z1 = reshape(data(:,9),nDistance,nSnr);
            x = 2:.2:10;
            y = 1:0.5:10;
            dfig;contourf(x,y,z1);
            title('% of 2 hits');
            dfig;contourf(reshape(data(:,10),nDistance,nSnr));
            
        end
        
        %% test Two Spot Two Kernel Off Start
        function testTwoSpotTwoKernelOffStart(obj)
            %n test on offStart and Signal-to-Noise ratio.
            nTest = 10^2;
            
            %real coordinate
            dataIdx = 1;
            nSnr = size(4:.2:10,2);
            data = NaN(size(nSnr,2),5);
            fullMatrixMean = zeros(1,2);
            for snr = 4:.2:10
                fprintf('SNR : %.3f\n',snr);
                
                %coordList ->output coord
                coordList = NaN(nTest,4);
                %coordList -> starting coordinate
                inputCoordList = NaN(nTest,4);
                groundTruthCoordList = NaN(nTest,4);
                goneCoord = NaN(nTest,4);
                nHit1 = 0;
                nHit2 = 0;
                for cTest = 1:nTest
                    %create test image
                    imgSize = [ 64 64];
                    coordinate = [ 32 32;32 38];
                    groundTruthCoordList(cTest,1:2) = coordinate(1,:);
                    groundTruthCoordList(cTest,3:4) = coordinate(2,:);
                    featureSize = 2;
                    [nseImg,~,~] = createMixtureFitTestImage('coordinates',coordinate,...
                        'snr',snr,'imgSize',imgSize,'featureSize',featureSize);
                    offSet = (rand(2,2)-[0.5,0.5;0.5,0.5])*8;
                    inputCoordList(cTest,1:2) = (coordinate(1,:)+offSet(1,:));
                    inputCoordList(cTest,3:4) = (coordinate(2,:)+offSet(2,:));
                    [newCoord,~,~,fullMatrix] = spotMMFit(nseImg,(coordinate+offSet),featureSize,'fitNPlusOne',1,'amplitudeCutOff',0.001);
                    if ~(isempty(newCoord))
                        if size(newCoord,1) == 2
                            coordList(cTest,1:2) = newCoord(1,:);
                            coordList(cTest,3:4) = newCoord(2,:);
                            nHit2 = nHit2 + 1;
                        elseif size(newCoord,1) == 1
                            coordList(cTest,1:2) = newCoord(1,:);
                            nHit1 = nHit1 + 1;
                        end
                        
                    else
                        goneCoord(cTest,:) = inputCoordList(cTest,:);
                        coordList(cTest,:) = NaN;
                    end
                end
                
                dfig,plot(coordList(:,1),coordList(:,2),'ob');
                title((snr));
                hold on
                plot(groundTruthCoordList(:,1),groundTruthCoordList(:,2),'xg');
                plot(groundTruthCoordList(:,3),groundTruthCoordList(:,4),'xg');
                
                plot(inputCoordList(:,1),inputCoordList(:,2),'or');
                plot(coordList(:,3),coordList(:,4),'ok');
                
                title((snr));
                plot(inputCoordList(:,3),inputCoordList(:,4),'or');
                axis equal;
                fprintf('Time found 1 coord : %d\nTime found 2 coord : %d\n',nHit1,nHit2);
            end
        end
        
        %% test N spot M kernel
        function testNSpotMKernel(obj)
            
            nSnr = size((8:.2:10),2);
            nSpot = 9;
            nTest = 10^2;
            nKernel = size((1:9),2);
            data = NaN(nSnr+nKernel,11);
            dataIdx = 1;
            for snr = 8:.2:10
                fprintf('snr : %.2f\n',snr);
                for m = 1:9
                    fprintf('Number of Kernel : %d\n',m);
                    coordList = NaN(nTest,18);
                    groundTruthCoordList = NaN(nTest,18);
                    nHit = zeros(1,9);
                    for cTest = 1:nTest
                        [nseImg,~,~,~,coordinate] = createNoiseImage('snr',snr,'cluster',1,'nCluster',1,'randomCoord',0);
                        [newCoord,~,~,fullMatrix] = spotMMFit(nseImg,coordinate(1:m,:),2,'fitNPlusOne',1);
                        
                        if size(newCoord,1) > 9
                            disp('AHHH');
                        else
                            nHit(1,size(newCoord,1)) = nHit(1,size(newCoord,1)) + 1;
                        end
                        for i = 1:size(newCoord,1)
                            coordList(cTest,(i-1)*2+1:(i-1)*2+2) = newCoord(i,:);
                        end
                        for i = 1:size(coordinate,1)
                            groundTruthCoordList(cTest,(i-1)*2+1:(i-1)*2+2) = coordinate(i,:);
                        end
                    end
                    dfig;
                    hold on;
                    for row = 1:size(groundTruthCoordList,2)/2
                        plot(groundTruthCoordList(1,(row-1)*2+1),groundTruthCoordList(1,(row-1)*2+2),'xr');
                    end
                    for row = 1:size(coordList,2)/2
                        plot(coordList(:,(row-1)*2+1),coordList(:,(row-1)*2+2),'ob');
                    end
                    title(m);
                    
                    data(dataIdx,1) = snr;
                    data(dataIdx,2) = m;
                    data(dataIdx,3:11) = nHit;
                    dataIdx = dataIdx + 1;
                end
            end
            
        end
        
        
        
        function testFit(obj)
            testIndex = 201;
            nKernel = 3;
            nSpot = 2;
            nSnr = size((2:.2:10),2);
            nTest = 10^2;
            imgSize = 128;
            featureSize = 3;
            
            coordinate = [ [1 1] * imgSize/2 - featureSize*0.7*1.3;...
                [1 1] * imgSize/2 + featureSize * 0.7*1.3];
            %coordinate = [ [1 1] * imgSize/2 - featureSize*0.7*1.3;...
            %   [1 1] * imgSize/2 + featureSize * 0.7*1.3; [1 1] * imgSize/2 + 3 * featureSize * 0.7 * 1.3];
            goodnessOfFit = [];
            for cKernel = 1:nKernel
                fprintf('nKernel : %d\n',cKernel);
                for snr = 2:1:10
                    fprintf('snr : %.1f\n',snr);
                    for cTest = 1:nTest
                        fprintf('test #%d\n',cTest);
                        [nseImg,~,~] = createMixtureFitTestImage('coordinates',coordinate,...
                            'snr',snr,'imgSize',[imgSize imgSize],'featureSize',featureSize);
                        if cKernel < 3
                            [~,~,~,~,~,stats] = spotMMFit(nseImg,coordinate(1:cKernel,:),2,'debug',1,'fitNPlusOne',1);
                        else
                            inCoord = coordinate;
                            inCoord(end+1,:) = coordinate(end,:) + rand(1,2)*5;
                            [~,~,~,~,~,stats] = spotMMFit(nseImg,inCoord,2,'debug',1,'fitNPlusOne',1);
                        end
                        if ~isempty(stats)
                            st = cat(2,stats,zeros(size(stats,1),3));
                            st(:,7) = testIndex;
                            st(:,8) = snr;
                            nFitSpot = size(stats,1);
                            if cKernel == 1
                                for cSpot = 1:nFitSpot
                                    if cSpot < 2
                                        st(cSpot,6) = 1;
                                    end
                                end
                                %                           elseif cKernel == 2
                                %                               for cSpot = 1:nFitSpot
                                %                                   if cSpot < 2
                                %                                       st(cSpot,6) = 1;
                                %                                   end
                                %                               end
                            end
                        else
                            st = [];
                        end
                        goodnessOfFit = cat(1,goodnessOfFit,st);
                    end
                end
                testIndex = testIndex + 1;
            end
            save('goodnessOfFit','goodnessOfFit');
        end
        
    end %methods
end % classdef% function out = setup
