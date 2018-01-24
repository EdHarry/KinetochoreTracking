function [adepth,algorithm,res,msg]=readTiffOptions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  READS TIFF OPTIONS FROM UIVIEWPANEL 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Opens/activates ViewPanel
auxV = findobj(0,'Type','figure','Tag','UIVIEWPANEL');
if isempty(auxV)
   uiviewpanel;
   % Now look for this object
   auxV = findobj(0,'Type','figure','Tag','UIVIEWPANEL');
else
   figure(auxV);
end

% Reads current Tiff Options
ld=findobj(auxV,'Tag','UIVIEWMENU_FILE_OPEN');
adal=get(ld,'UserData');

% Restoring adepth information
switch adal(1);
case '1'
   adepth=-1;	% One of the search algorthm selected
case '2'
   adepth=8;
case '3'
   adepth=9;
case '4'
   adepth=10;
case '5'
   adepth=11;
case '6'
   adepth=12;
case '7'
   adepth=13;
case '8'
   adepth=14;
case '9'
   adepth=15;
case 'A'
   adepth=16;
end

% Restoring algorithm information
switch adal(2)
case '0'
   algorithm=-1;	% No search algorithm   
case '1'
   algorithm=0;	% Simple algorithm
case '2'
   algorithm=1;	% Smart algorithm
end

% Before calling imreadnd, the parameters 'res' and 'msg' must receive a default value, in case no file is chosen for loading 
res=2;	%If res stays equal to 2, that means that no file has been chosen for load (imreadnd returns either 0 or -1)
msg='';
