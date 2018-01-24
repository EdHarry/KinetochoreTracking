function [firstFileName]=callGetStackBoundaries

% Select first name
[fName,dirName] = uigetfile('*.tif','View Panel: Open ...');
if(isa(fName,'char') & isa(dirName,'char'))
   firstfilename=strcat(dirName,fName);
   cd(dirName);
end;

% Resets string
displ('Evaluating... please wait');

% Looks for global min and max intensity values over the stack
[gmin,gmax,fileMin,fileMax]=getStackBoundaries(firstfilename);

% Calculates boundaries in bit
if gmin==0
   bmin=0;
else
   bmin=fix(log2(gmin));
end
if gmax~=0
   if (log2(gmax)-fix(log2(gmax)))>0
      bmax=fix(log2(gmax))+1;
   else
      bmax=log2(gmax);   
   end
else
   xmax=0;
end

% Creates output string
msg=strcat('Result: global min= ',num2str(gmin),' (->',num2str(bmin),'bit)',' / global max= ',num2str(gmax),' (->',num2str(bmax),'bit)');
%msg=strcat('Result: global min= ',num2str(gmin),'(->',num2str(bmin),'bit)','[',fileMin,'] / global max= ',num2str(gmax),'(->',num2str(bmax),'bit)',' [',fileMax,']');

% Displays string in message window
displ(msg);








