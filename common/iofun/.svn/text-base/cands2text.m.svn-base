 defaultFileName = '*.mat';    
 [firstCandsName, firstCandsPath] = uigetfile('.mat', 'Please select the first speckle/vesicle data file', defaultFileName);

 
 
 firstCandsName = [firstCandsPath, firstCandsName];

 candsStackList = getFileStackNames(firstCandsName);
 
 cd(firstCandsPath);
 
 allCoordinates = [];
 for i = 1 : length(candsStackList)
     tempFileName = char(candsStackList(i));
     temp = load(tempFileName);
     allCoordinates = [];

    if (i < 10)
         tempFileName = ['coordinates', '000', num2str(i), '.txt'];
     elseif (i < 100)
         tempFileName = ['coordinates', '00', num2str(i), '.txt']; 
     elseif (i < 1000)
         tempFileName = ['coordinates', '0', num2str(i), '.txt'];
     else
         tempFileName = ['coordinates', num2str(i), '.txt'];
    end
    
    tempFileName = [firstCandsPath, tempFileName];
     fid = fopen(tempFileName, 'wt');
     
     for j = 1 : length(temp.candsSP)
         if (temp.candsSP(j).status == 1)
             fprintf(fid, '%10.3f  %10.3f\n', temp.candsSP(j).Lmax(2), temp.candsSP(j).Lmax(1));
         end
     end
     fclose(fid);
     
 end
 
 
