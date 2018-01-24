function [gamma,gammaO] = makiInlierSisterAutoCorrPole
% EHarry Jan 2012

c=0;
cO=0;

list = searchFiles('-makiData-',[],uigetdir);

for i = 1:length(list)
    dataStruct = makiLoadDataFile('MCAINSH',fullfile(list{i,2},list{i,1}));
    
    if isempty(dataStruct.sisterList_pole_aligned)
        continue
    end
    
    maxLag = dataStruct.dataProperties.movieSize(4)./8;
    unalignedAndLaggingSisters = [dataStruct.updatedClass_pole(1).sistersLagging dataStruct.updatedClass_pole(1).sistersUnaligned];
    unalignedAndLaggingSistersO = [dataStruct.updatedClass(1).sistersLagging dataStruct.updatedClass(1).sistersUnaligned];
    sisterList = dataStruct.sisterList_pole_aligned;
    sisterListO = dataStruct.sisterList;
    for iSis = 1:length(sisterList)
        if ~ismember(iSis,unalignedAndLaggingSisters)
            if ~isempty(sisterList(iSis).distances)
                c=c+1;
                centre = (sisterList(iSis).coords1(1:4:end,1) + sisterList(iSis).coords2(1:4:end,1))./2;
                centreStd = 0.5.*sqrt((sisterList(iSis).coords1(1:4:end,4)).^2 + (sisterList(iSis).coords2(1:4:end,4)).^2);
                centreDiff(:,1) = diff(centre);
                centreDiff(:,2) = sqrt( sum( [centreStd(2:end) centreStd(1:end-1)].^2 ,2) );
                traj(c).observations = centreDiff;
                clear centreDiff
            end
        end
    end
    
    for iSis = 1:length(sisterListO)
        if ~ismember(iSis,unalignedAndLaggingSistersO)
            if ~isempty(sisterListO(iSis).distances)
                cO=cO+1;
                centre = (sisterListO(iSis).coords1(1:4:end,1) + sisterListO(iSis).coords2(1:4:end,1))./2;
                centreStd = 0.5.*sqrt((sisterListO(iSis).coords1(1:4:end,4)).^2 + (sisterListO(iSis).coords2(1:4:end,4)).^2);
                centreDiff(:,1) = diff(centre);
                centreDiff(:,2) = sqrt( sum( [centreStd(2:end) centreStd(1:end-1)].^2 ,2) );
                trajO(cO).observations = centreDiff;
                clear centreDiff
            end
        end
    end
end

gamma = autoCorr(traj,maxLag);
gammaO = autoCorr(trajO,maxLag);

end

