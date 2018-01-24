function clrmp=wrappedColors

% wrappedColors is a wrap-around colormap 
% BLUE CYAN GREEN YELLOW RED MAGENTA BLUE
% intended to be used particularly for Orientation Data
% 
% SYNOPSIS   clrmp=ourColors 
% 
% INPUT      none     
%
% OUTPUT     a wrap-around colormap 
%
% DEPENDENCES   wrappedColors uses {  }
%               wrappedColors is used by { textureFilter }

% Alexnadre Matov, February 20th, 2003
 
R=zeros(180,1);
R(61:90,1)=[0:1/29:1]';
R(91:150,1)=1;
R(151:180,1)=[1:-1/29:0]';

G=zeros(180,1);
G(1:30,1)=[0:1/29:1]';
G(31:90,1)=1;
G(91:120,1)=[1:-1/29:0]';

B=zeros(180,1);
B(1:30,1)=1;
B(31:60,1)=[1:-1/29:0]';
B(121:150,1)=[0:1/29:1]';
B(151:180,1)=1;

clrmp=[R G B];
clrmp=[[0 0 0];clrmp];