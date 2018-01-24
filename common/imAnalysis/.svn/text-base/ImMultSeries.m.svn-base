function multImgSeries(f)
BIT_DEPTH = 2^16;

firstfilename = 'S:\scripps\Staurosporin\stauroRatio\ratio01.tif'
[filelist]=getFileStackNames(firstfilename);
fileName=char(filelist(1));
[fpath,fname] = fileparts(fileName);

mkdir([fpath filesep 'CorrectedImg']);

for i=1:length(filelist)
    
    
    fileName=char(filelist(i));
    %img_org=imreadnd2(fileName,0,BIT_DEPTH);
    img_org=imread(fileName,'tif');
    img_cor = uint16(f(i).* double(img_org));
    
    
    [fpath,fname] = fileparts(fileName);
    imwrite(img_cor,[fpath filesep 'CorrectedImg'  filesep 'corr_' fname '.tif'],'tif');
end