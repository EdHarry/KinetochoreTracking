function [pass,fail] = testFtestVSs2n
%TESTFTESTVSS2N Summary of this function goes here
%   Detailed explanation goes here


pass = [];
fail = [];
count = 0;

psf = [1 1 1];


for k = 1:5
    
    spot = [6 6 6 1];
    
    for i = 1:10
        nSpot = size(spot,1);
        spotNew = spot(end,:);
        spotNew = spotNew + [3 0 0 0];
        spot = [spot; spotNew];
        [~,imageO] = placeGausianSpots( spot , [psf(1) psf(3)], [11+5*(nSpot-1) 11 11]);
        for j = 0:100
            image = imageO + j*0.003*randn(size(imageO));
            [~,~,~,~,~,stats] = spotMMFit(image,spot(:,1:3),psf,'fitNPlusOne',1,'debug',1,'F_test_prob',0.05);
            count = count + 1;
            if ~isempty(stats)
                %pass = [pass; nSpot stats(1,6) stats(1,7)];
                %if size(stats,1) > 1
                for m = 1:size(stats,1)
                    fail = [fail; nSpot stats(m,6) stats(m,7)];
                end
                waitbar(count/(10*101*5));
            end
        end
    end
end

