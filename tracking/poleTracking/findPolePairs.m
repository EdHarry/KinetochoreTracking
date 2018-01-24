function dataStruct = findPolePairs( dataStruct )

if ~isfield(dataStruct,'sisterList') || isempty(dataStruct.sisterList)
    return
end

nFrames = dataStruct.dataProperties.movieSize(4);

previousPoles = struct('pole1',zeros(2,nFrames),'pole2',zeros(2,nFrames));

loop = true;
while loop
    dataStruct = makiFindAllPoles(dataStruct);
    
    nTracks = length(dataStruct.tracks);
    nSis = length(dataStruct.sisterList);
    sList = dataStruct.sisterList(1).trackPairs;
    sListColumns = size(sList,2);
    
    poles = dataStruct.poles;
    
    if isempty(poles) || (length(poles.pole1Track) < 2 && length(poles.pole2Track) < 2)
        loop = false;
    else
        repeatedPoles = false;
        counter = 0;
        while ~repeatedPoles && counter < length(previousPoles)
            pole1Found = false;
            pole2Found = false;
            counter = counter + 1;
            if length(poles.pole1Track) == size(previousPoles(counter).pole1,1) && length(poles.pole2Track) == size(previousPoles(counter).pole2,1)
                if length(poles.pole1Track) == 1 && all(getFeatIdxs(poles.pole1Track) == previousPoles(counter).pole1)
                    pole1Found = true;
                elseif length(poles.pole1Track) == 2 && ((all(getFeatIdxs(poles.pole1Track(1)) == previousPoles(counter).pole1(1,:)) && all(getFeatIdxs(poles.pole1Track(2)) == previousPoles(counter).pole1(2,:))) || (all(getFeatIdxs(poles.pole1Track(2)) == previousPoles(counter).pole1(1,:)) && all(getFeatIdxs(poles.pole1Track(1)) == previousPoles(counter).pole1(2,:))))
                    pole1Found = true;
                end
                
                if length(poles.pole2Track) == 1 && all(getFeatIdxs(poles.pole2Track) == previousPoles(counter).pole2)
                    pole2Found = true;
                elseif length(poles.pole2Track) == 2 && ((all(getFeatIdxs(poles.pole2Track(1)) == previousPoles(counter).pole2(1,:)) && all(getFeatIdxs(poles.pole2Track(2)) == previousPoles(counter).pole2(2,:))) || (all(getFeatIdxs(poles.pole2Track(2)) == previousPoles(counter).pole2(1,:)) && all(getFeatIdxs(poles.pole2Track(1)) == previousPoles(counter).pole2(2,:))))
                    pole2Found = true;
                end
            end
            if pole1Found && pole2Found
                repeatedPoles = true;
            end
        end
        
        if ~repeatedPoles
            if length(poles.pole1Track) == 1
                pole1Tmp = getFeatIdxs(poles.pole1Track);
            else
                pole1Tmp = [getFeatIdxs(poles.pole1Track(1)) ; getFeatIdxs(poles.pole1Track(2))];
            end
            
            if length(poles.pole2Track) == 1
                pole2Tmp = getFeatIdxs(poles.pole2Track);
            else
                pole2Tmp = [getFeatIdxs(poles.pole2Track(1)) ; getFeatIdxs(poles.pole2Track(2))];
            end
            
            previousPoles = [previousPoles struct('pole1',pole1Tmp,'pole2',pole2Tmp)];%#ok<AGROW>
            
            if length(poles.pole1Track) == 2
                sisIdx = zeros(1,sListColumns);
                for i = 1:2
                    if poles.pole1TrackIdx(i) > 0
                        sisIdx(i) = poles.pole1TrackIdx(i);
                    else
                        dataStruct.tracks(nTracks+1) = poles.pole1Track(i);
                        nTracks = nTracks + 1;
                        sisIdx(i) = nTracks;
                    end
                end
                dataStruct.sisterList(nSis+1).distances = 0;
                nSis = nSis + 1;
                sList = [sList; sisIdx];%#ok<AGROW>
            end
            
            if length(poles.pole2Track) == 2
                sisIdx = zeros(1,sListColumns);
                for i = 1:2
                    if poles.pole2TrackIdx(i) > 0
                        sisIdx(i) = poles.pole2TrackIdx(i);
                    else
                        dataStruct.tracks(nTracks+1) = poles.pole2Track(i);
                        nTracks = nTracks + 1;
                        sisIdx(i) = nTracks;
                    end
                end
                dataStruct.sisterList(nSis+1).distances = 0;
                sList = [sList; sisIdx];%#ok<AGROW>
            end
            
            dataStruct.sisterList(1).trackPairs = sList;
            
            dataStruct = KTPairTracking(dataStruct);
        else
            loop = false;
        end
    end
end


    function featIdx = getFeatIdxs(track)
        featIdx = zeros(1,nFrames);
        startTime = track.seqOfEvents(1,1);
        endTime = track.seqOfEvents(2,1);
        featIdx(startTime:endTime) = track.tracksFeatIndxCG;
    end

end

