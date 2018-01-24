function plotspec(xresolfinal, xresolintermediate) 
%PLOTSPEC 	service 3DPlot data from fluorimeter.
%
% SYNOPSIS	plotspec(xresolfinal, xresolintermediate)
%	
% INPUT 		xresolfinal : resolution in x of the graph
%				xresolintermediate : resolution in x of the interpolation

%				default value : xresolfinal = lowest step between data sets
%				If xresolfinal is set, xresolintermediate is 1/3 xresolfinal
%				otherwise, it is 1/3 lowest step between data sets

% Initialization of variables 
	ztemp=[];
	zmax=1;
	oldpath = pwd;
	ycoor=[];
	zcoor=[];
	newzcoor=[];
	em=[];
   wavelength=[];

% Selection of the first and last files 
	[firstx,filepath]=uigetfile('*.prn', 'Select First Spectometer File');
	cd(filepath)
	[lastx,filepath]=uigetfile('*.prn', 'Select Last Spectometer File');

% Assess the filname, the starting and ending excitation wave length  
	[path,name,ext,dummy] = fileparts(firstx);
   startx=str2num(name((size(name,2)-2):(size(name,2))));
   fname=name(1:(size(name,2)-3));
	[path,name,ext,dummy] = fileparts(lastx);
	endx=str2num(name((size(name,2)-2):(size(name,2))));

% Assess the list and the number of files to be read
	dirlist=dir(strcat(fname,'*',ext));
   nbfiles=length(dirlist);
   
% Reading the data and filling the matrix
	j=1;
	oldx=0;
	minxstep=endx;
	for dlnb=1:nbfiles
   	[path,name,ext,dummy] = fileparts(dirlist(dlnb).name);
   	x=str2num(name((size(name,2)-2):(size(name,2))));
   	if ((x>=startx)&(x<=endx))&(size(name(1:(size(name,2)-3)))==size(fname))
         if name(1:(size(name,2)-3))==fname
            if x-oldx<minxstep 							%check for the minimum resolution in x
         		minxstep=x-oldx;
      		end
      		a=load(strcat(fname,num2str(x),ext));  % Data are stored in a temporary matrix a[2,nb of data]
         	ztemp=cat(1,ztemp,a(:,2));  				% The intensities are stored and catenated in the ztemp matrix [1,inf]
         	em(j,1)=x;										% The matrix em is filled (excitationwl, emmissionwl(start),_
         	em(j,2)=a(1,1);								% emmissionwl(end), step, numb step)
   			em(j,3)=a(size(a,1),1);
   			em(j,4)=(a(2,1)-a(1,1));					% Ystep is assessed in the matrix
   			em(j,5)=(em(j,3)-em(j,2))/em(j,4)+1;	%
      		j=j+1;
         	oldx=x;
         end
   	end
   end
   
% Interpolation of the data (z) for the biggest nb of data into the zcoor matrix
	j=0;
   for i=1:size(em,1) 
      zcoor=cat(1,zcoor,interp1((em(i,2):em(i,4):em(i,3)),(ztemp(j+1:j+em(i,5))'),(em(i,2):((em(i,3)-em(i,2))/(max(em(:,5))-1)):em(i,3)),'cubic'));
      j=j+em(i,5);
   end
   
% Assess what parameters have been input
	if nargin==0
   	xresolfinal=minxstep;
   	xresolintermediate=(minxstep/3);
	end
	if nargin==1
   	xresolintermediate=(xresolfinal/3);
	end
   
% First interpolation - Xresolution
	xcoor=startx:xresolintermediate:endx;
	newymin=interp1(em(:,1), em(:,2), xcoor);
   newymax=interp1(em(:,1), em(:,3), xcoor);
   zcoor=interp2(em(:,1),(0:(1/(max(em(:,5))-1)):1),zcoor',xcoor,(0:(1/(max(em(:,5))-1)):1)','cubic');
   maxnewymax=(max(newymax));
   minnewymin=(min(newymin));

% Formation of the final z matrix
	j=1;
   for i= startx:xresolintermediate:endx
      ytemp=newymin(j):((newymax(j)-newymin(j))/(max(em(:,5))-1)):newymax(j);
      newytemp=minnewymin:((maxnewymax-minnewymin)/(max(em(:,5))-1)):maxnewymax;
      newzcoor=cat(1,newzcoor,interp1(ytemp,zcoor(:,j),newytemp));
      j=j+1;
   end 
   for a=1:size(newzcoor,1)
      for b=1:size(newzcoor,2)
         if isnan(newzcoor(a,b))==1
            newzcoor(a,b)=0;
         end
      end
   end
   newzcoor=newzcoor';
   
% Final interpolation for X given resolution
	newzcoor=interp2(xcoor,newytemp',newzcoor,startx:xresolfinal:endx,newytemp','nearest');
   
% Plot the surface
	zmax=max(max(newzcoor));
	axis([startx endx minnewymin maxnewymax 0 zmax]);
	colormap(jet);
   surface(startx:xresolfinal:endx, newytemp, newzcoor);
   shading interp; 

% change back the path to old path
	cd(oldpath);