function [G,W,S,mnWindows]=framework(i,g,w)
% FRAMEWORK creates a 2D regular grid
%
% SYNOPSIS  [G,W,S,nWindows]=framework(i,g,w)
%
% INPUT   i : image size  [iy ix]
%         g : grid size   [gy gx]
%             dimensions must be odd
%         w : window size [wy wx] (optional)
%             (dimensions must be odd - THIS CHECK HAS BEEN DISABLED - see lines
%             35 - 39)             
%             the window size must be equal or smaller than the grid size
%
% OUTPUT  G : matrix with all grid centers [Gy Gx]
%         W : matrix with all region edge coordinates [Wy0 Wx0 Wy Wx] for correlation template
%         S : matrix with all region edge coordinates [Sy0 Sx0 Sy Sx] for source template 
% mnWindows : number of subwindows in y and x direction
%
% DEPENDENCES      framework uses { }
%                  framework is used by { selectRegion }
%
% Aaron Ponti, 2001

% Check input parameters
if nargin<2
   error('Not enough input arguments defined');
elseif nargin==2
   w=g;
elseif nargin==3
   if isempty(w)
      w=g;        
   end
end

% Masks dimensions have to be odd
% fg=find(~mod(g,2));
% g(fg)=g(fg)+1;
% fw=find(~mod(w,2));
% w(fw)=w(fw)+1;

% Image must be larger than both masks
y=max(w(1),g(1));
x=max(w(2),g(2));
if i(1)<y | i(2)<x
   error('One of the masks is too large');
end

% Calculate first center
g0=[fix(g(1)/2)+1 fix(g(2)/2)+1];   

% Number of grid centers in y and x direction
nY=1:fix(i(1)/g(1));
nX=1:fix(i(2)/g(2));

% Their coordinates
cy=(g0(1)+(nY-1)*g(1));
cx=(g0(2)+(nX-1)*g(2));

% Vertical coordinates
for j=1:length(nX)
   Gx(1:length(nY),j)=cx(j);
end
Gx=Gx(:);
for k=1:length(nX)
   Gy((k-1)*(length(nY))+1:k*(length(nY)))=cy;
end
Gy=Gy(:);
% Matrix G(Gy Gx)
G=Gy;
G(1:length(Gx),2)=Gx;

% Minimum and maximum allowed grid centers
minY=fix(w(1)/2)+1;
maxY=i(1)-fix(w(1)/2);
minX=fix(w(2)/2)+1;
maxX=i(2)-fix(w(2)/2);

% Number of windows in y and x direction
p=find(cy<minY | cy>(maxY+1));
if ~isempty(p)
    numberY=length(nY)-length(p);
else
    numberY=length(nY);
end
clear p;

q=find(cx<minX | cx>(maxX+1));
if ~isempty(q)
    numberX=length(nX)-length(q);
else
    numberX=length(nX);
end
clear p;
 
mnWindows=[numberY numberX];

% Remove all grid centers from G which are to close to the borders
p=find(G(:,1)<minY | G(:,1)>(maxY+1));
if ~isempty(p)
   G(p,:)=[];
end
q=find(G(:,2)<minX | G(:,2)>(maxX+1));
if ~isempty(q)
   G(q,:)=[];
end

% TEMP
G=sortrows(G,1);

% Calculate region vertex coordinates (correlation template W and source template S)
W(:,1)=G(:,1)-fix(w(1)/2);
W(:,2)=G(:,2)-fix(w(2)/2);
S(:,1)=G(:,1)-fix(w(1)/4);
S(:,2)=G(:,2)-fix(w(2)/4);

if mod(w(1),2)==0
    W(:,3)=G(:,1)+fix(w(1)/2)-1; % To avoid overlap with the first row of the next window
    S(:,3)=G(:,1)+fix(w(1)/4)-1;
else
    W(:,3)=G(:,1)+fix(w(1)/2); % No overlap anyway
    S(:,3)=G(:,1)+fix(w(1)/4)-1;
end
if mod(w(2),2)==0
    W(:,4)=G(:,2)+fix(w(2)/2)-1;
    S(:,4)=G(:,2)+fix(w(2)/4)-1;  
else
    W(:,4)=G(:,2)+fix(w(2)/2); 
    S(:,4)=G(:,2)+fix(w(2)/4);
end


