function uiSSDTrackPanelWriteLog(indx,cde,pos,posRef,shape,lori,nPix,fname)
% service function for the SSD Track panel appending the 
% the results of a match to the log file

fid = fopen(fname,'a');
if(fid <0)
   return;
end;

if(~isempty(indx))
   fprintf(fid,'%6d',indx);
end
fprintf(fid,'%6d',cde);
fprintf(fid,'%10.5f %10.5f',pos+(shape*posRef')');
fprintf(fid,'%10.7f %10.7f %10.7f %10.7f ',shape');

% make SVD of the shape matrix
% is not of much interest right now GD Nov-09-1998
%[u,s,v] = svd(shape);
%fprintf(fid,'%10.7f %10.7f %10.7f %10.7f',...
%   atan2(u(3),u(1)),atan2(v(3),v(1)),s(1),s(4));
if(isempty(lori))
   lori = 0.0;
end;
fprintf(fid,'%10.7f ',lori);

fprintf(fid,'%6d',nPix);

fprintf(fid,'\n');
fclose(fid);


