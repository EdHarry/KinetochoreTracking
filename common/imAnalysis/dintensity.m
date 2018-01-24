function [aint,dint,rint]=dintensity(stack,displ,region,noi)
% DINTENSITY sums the intensities of all pixels in a defined region for all images in a stack
%  and compares them. 
%
%  SYNOPSIS [aint,dint,rint]=DINTENSITY(stack,displ,region,noi)
%
%  INPUT    stack     Stack created using the IMREADSTACK function
%           displ     Display of results. Options are: 0 - none
%                                                      1 - text
%                                                      2 - graphical
%                                                      3 - text and graphical
%           region    (Optional) Matrix (1,4) defining the region of interest [x0 y0 x y]
%                                x must be >= x0 and y must be >=y0
%                                If region is not defined the entire image will be considered
%           noi       (Optional) Number of images from the stack to be analyzed
%                                If noi is not defined all images in the stack will be considered
%
%  OUTPUT   aint      Sum of all intensities within the regions of interest
%           dint      Differences between overall intensities of one image in the stack 
%                     compared to the previous one
%           rint      Relative difference between two subsequent frames
%
%                     These vectors can be displayed in a table form or plotted

% *************************
%
% INPUT PARAMETERS CHECK
%
% *************************

% Check for input parameters existance
if nargin<1
   error('the input parameter stack is not optional!');
end
if nargin==1
   displ=0;
   region=[1 1 size(stack,1) size(stack,2)];
   noi=size(stack,3);
end
if nargin==2
   region=[1 1 size(stack,1) size(stack,2)];
   noi=size(stack,3);
end
if nargin==3
   noi=size(stack,3);
end

% Check for stack type (must be double for operations to be possible) 
if class(stack)~='double'
   stack=double(stack);
end
% Check for region of interest consistency (first)
if ~(size(region,1)==1 & size(region,2)==4) 
	error('region does not define a rectangular area');   
end
% Check for region of interest consistency (second)
if (region(1)<1) | (region(2)<1) | (region(3)>size(stack,1)) | (region(4)>size(stack,2))
   error('the selected region is not valid');   
end
% Check for region of interest consistency (third)
if ((region(3)<region(1)) | region(4)<region(2))
   error('region coordinates [x0 y0 x y] are not valid');
end
% Check for number of images
if noi>size(stack,3)
   error('stack does not contain so many images');
end
% Check for display option
if displ<0 | displ>4
   displ=0;
end

% *************************
%
% EVALUATION
%
% *************************

% Creating a vector 'aint' to store all total intensities of all regions from the stack
aint=zeros(noi,1);
% Creating a vector 'dint' to store all intensity differences between the regions of 
% interest of to sequential images from the stack
dint=zeros(noi-1,1);

% Calculating all total intensities from all regions of the stack
for counter1=1:noi
   % Initializing the current global intensity
   gint=0;
   % Calculating the global intensity over the entire region
   for counter2=region(1,2):region(1,4)
      for counter3=region(1,1):region(1,3)
         gint=gint+stack(counter2,counter3,counter1);
      end
   end
   % Storing the total intensity of the current region
   aint(counter1,1)=gint;
   % Storing the difference in intensity compared to the last image
   if counter1>1
      dint(counter1-1,1)=gint-aint(counter1-1,1);   
   end
end

% Calculating vector of relative difference intensities rint
rint=dint(1:noi-1)./aint(2:noi);

% *************************
%
% RESULTS OUTPUT
%
% *************************

% default switches
tswitch=0;
gswitch=0;

% display option chosen
switch displ
case 0
case 1
   tswitch=1;
case 2
   gswitch=1;
case 3
   tswitch=1;
   gswitch=1;
otherwise
end

% Display text
if tswitch==1
   fprintf(1,'\n\nOUTPUT\n------\n\n');
   fprintf(1,'   #  totInt   dInt    rInt\n---------------------------------\n');              
   for counter1=1:size(aint,1)
      if counter1<2
         fprintf(1,' %3.0d: %6.3f\n',counter1,aint(counter1,1));
      else
         fprintf(1,' %3.0d: %6.3f (%6.3f) %6.3f%%\n',counter1,aint(counter1,1),dint(counter1-1,1),dint(counter1-1,1)/aint(counter1)*100);
      end
   end
   fprintf(1,'\n');
end

% Display plots
if gswitch==1
   figure;
   H=plot(aint,'-s');
   title('Intensity vs Frame');
   figure;
   bar(rint*100);
   title('Relative intensity difference (%) vs dFrame');
   fprintf(1,'\n\nTwo plots generated.\n');
end


   
   



