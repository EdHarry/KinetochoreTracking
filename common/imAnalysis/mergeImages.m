function mergeImages

% mergeImages takes two sequences of images from two chosen by the user folders 
% and merges them horizontally
%
%
% SYNOPSIS   mergeImages
%
% INPUT      none        :   manual choise
%             
%
% OUTPUT     none        :   writes to disk
%
%
% DEPENDENCES   mergeImages uses {  }
%               mergeImages is used by {  }
%
% Alexandre Matov, January 7th, 2003

% choose first input directory and its first image
[file1Name,dir1Name] = uigetfile('*.tif','Choose first input directory and its first image');


if(isa(file1Name,'char') & isa(dir1Name,'char'))
   % Recover all file names from the stack
   outFileList1=defineStackNames([dir1Name file1Name]);
   % Number of files 
   n1=length(outFileList1);
else
   return
end

aux1=imfinfo([dir1Name,filesep,file1Name]);
BitDepth1=aux1.BitDepth;
if aux1.ColorType~='grayscale'
    error('Not a Gray Scale Image');
end

[y1,x1]=size(imread([dir1Name,filesep,file1Name]));
 
% choose second input directory and its first image
[file2Name,dir2Name] = uigetfile('*.tif','Choose second input directory and its first image');

if(isa(file2Name,'char') & isa(dir2Name,'char'))
   % Recover all file names from the stack
   outFileList2=defineStackNames([dir2Name file2Name]);
   % Number of files 
   n2=length(outFileList2);
else
   return
end

aux2=imfinfo([dir2Name,filesep,file2Name]);
BitDepth2=aux2.BitDepth;
if aux2.ColorType~='grayscale'
    error('Not a Gray Scale Image');
end

[y2,x2]=size(imread([dir2Name,filesep,file2Name]));


y=min(y1,y2);
x=x1+x2;
BitDepth=max(BitDepth1,BitDepth2);
n=min(n1,n2);
s=length(num2str(n));
strg=sprintf('%%.%dd',s); 

% choose output directory and file name
[file3Name,dir3Name] = uiputfile('Specify a file name','Choose output directory and file name');


h=waitbar(0,'Please wait! The program is merging the images');
for i=1:n

    % Crruent index
    indxStr=sprintf(strg,i);
    
    % Current filenames
    currentFile1=char(outFileList1(i));
    currentFile2=char(outFileList2(i));
   
    % Read images from disk
    I1=imread(currentFile1);
    I2=imread(currentFile2);
    
    % Merge images
    Ir=zeros(y,x);
    Ir(1:y,1:x1)=I1;
    Ir(1:y,(x1+1):x)=I2;
    
    % Convert image to integer
    switch BitDepth
    case 8
        Ir=uint8(Ir);
    case 16
        Ir=uint16(Ir);
    otherwise
        error('BitDepth Not Supported');
    end
    
    % Write resulting image to disk
    imwrite(Ir,[dir3Name,filesep,file3Name,indxStr,'.tif']);
    
    waitbar(i/n,h);  
end
close(h);


