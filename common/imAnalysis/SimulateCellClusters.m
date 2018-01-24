function[mpm]=SimulateCellClusters(m1,cbeh,cspid,np,dt,matsiz,celldia)
% SimClusterMatVary simulates a time series (yielding an mpm file) of 
% moving cells; it can simulate different kinds of behaviour, including
% cell division, constant or varying cell speed, constant or varying degrees
% of cell attraction/adhesion, and varying cell size
%
% SYNOPSIS [mpm]=SimulateCellClusters(m1,cbeh,cspid,np,dt,matsiz,
%                   celldia)
%       
% INPUT      m1: initial distribution of objects (clustered or
%                non-clustered matrix), nx2 matrix containing the (x,y) 
%                coordinates of n points; m1 can be the first entry in
%                another mpm, such as MPMexperimental(:,1:2)
%            
%           cbeh: cell behaviour - this parameter means the relative
%               adhesion/attraction between cells, values should be between 1 and 0
%               if this parameter is 1,
%               then the cells stick to each other, and during their
%               moveement, they are restricted to positions in the vicinity
%               of other cells; if this parameter is 0, then cells are able
%               to move freely
%           
%            cspid: cell average speed (pixels per frame)
%
%            np: number of planes; this parameter denotes how many times 
%                the simulation is performed and how many planes the final
%                mpm matrix will have

%            dt: doubling time of the population - this parameter denotes 
%                how fast the cells divide. Is in units of planes, not
%                in seconds; thus this parameter should not be
%                significantly smaller than the number of planes np, or
%                else there will be an explosion of cells
%                if dt=0, there is no growth
%
%            matsiz: matrix/image size [msx msy] in pixels, e.g. [1344 1024]
%
%            celldia:  cell diameter, this parameter is optional, if there 
%               is no entry, the cell diameter is taken from the nearest 
%               neighbor profile and kept to this value for the entirety of 
%               the movie
%               values of this parameter should denote unitless 
%               multiples of the cell diameter as determined above; e.g. if
%               celldia is a vector [1 1 1 1 .... 2 2 2 2], it means that
%               the cell double their original diameter over the course of
%               the movie; for values below 1, the cells shrink
%               
%
%               IPORTANT NOTE: 
%               cbeh, celldia, cspid : these parameters can either be 
%               single values, then this value will remain constant for the
%               entire movie - or they can be vectors. 
%               As an example: cspid is a vector, of the form 
%               cspid=[2 2 2 2 3 3 4 4 4 4].
%               For the simulation, the speed is 2 pixels/frame in
%               frames 1-4, 3 pixels/frame in frames 5-6, and 4 pixels/frame
%               in frames 7-10. The vector doesn't necessarily need to have
%               the same length as the movie (which is  = np).
%               If cspid is shorter than np, then the missing values
%               will be filled up as the last entry in the vector. In the
%               above example, this means that if np=20, then the speed in
%               frames 11-20 will be 4 pixels/frame, as in frame 10.
%                           
%
% OUTPUT     mpm: matrix containing all the coordinates over time
%             
%
% DEPENDENCES   SimClusterMat  uses {distanceMatrix, nearNeiDistProbHisto}
%               SimClusterMat  is used by { }
%
% Dinah Loerke, September 13th, 2004
%               December 09, 2005



% initialize temporary veriables
% mtemp contains all points in a given frame
% movvec contains last frame-to-frame movement vector
mtemp=[nonzeros(m1(:,1)), nonzeros(m1(:,2)) ];
movvec=mtemp;
movvec(:)=0;
[m,n]=size(mtemp);
mpm=mtemp;

plot(mtemp(:,1),mtemp(:,2),'r.');
pause(0.01);
m2=mtemp;


%calculate growth rate from doubling time
%note: growth rate is equivalent for probability of a single object
%dividing
if (dt==0)
    disp(['no growth, constant number of ',num2str(m),' cells per frame']);
    gr=0;
else
    gr=log(2)/dt;
    disp(['growth rate (%) per plane = ',num2str(100*gr)]);
end

lastmov=zeros(m,1);


%for first plane of the distribution, determine adhesion profile, i.e. the
%nearest neighbor probability distribution
[pd,pdfilter]=nearNeiDistProbHisto(m1);
%the average cell diamater is the position of the maximum in the smoothed
%probability distribution
csiz_max = min ( find (pdfilter == max(pdfilter)) ) ;

%the minimum constitutes the repulsion distance - (cells don't approach
%closer than this distance
csiz_min = min (find(pdfilter(1:csiz_max)>0) );


h = waitbar(0,'simulation');

%since cell don't move completely randomly - on the time scales we're
%considering they tend to move in the same direction for a bit, the vector
%lastAng is used to remember what angle/direction all the cells moved in in
%the last frame
lastAng=2*pi*rand(m,1)-pi;

%loop over np (number of frames)

for inp=1:np
    waitbar(inp/np);
    
    % number of cells in this frame (can be changed through cell division
    % in the last frame)
    [nx,ny]=size(mtemp);
    nden=nx/(matsiz(1)*matsiz(2));
    %calculate nearest neighbour distance in densest hexagonal packing for
    %this density (is needed for edge effects)
    thnn=sqrt(2/(sqrt(3)*nden));
    
    %disp(['cycle number = ',num2str(inp),'  number of points = ',num2str(nx)]);
    
        
    %=====================================================================
    % in each frame, determine applicable speed, adhesion potential profile 
    % and cell size
    %=====================================================================
    
    %SPEED
    % if there's only one value, assume constant cspid for all planes
    if(length(cspid)==1)
        spid=cspid;
    else
    % if there's more than one value, but the vector is not as long as the
    % number of planes reached, take the last value of cspid vector
        if (length(cspid)<inp)
           spid=cspid(length(cspid)); 
        else
           spid=cspid(inp); 
        end
    end

    %ADEHESION BEHAVIOR
    % if there's only one value, assume constant cbeh for all planes
    %cellbeh is zeros and ones to determine adhesion or lack thereof
    if(length(cbeh)==1)
        beh=cbeh;
    else
    % if there's more than one value, but the vector is not as long as the
    % number of planes reached, take the last value of cspid vector
        if (length(cbeh)<inp)
           beh=cbeh(length(cbeh)); 
        else
           beh=cbeh(inp); 
        end
    end

    %CELL SIZE
    % if there's only one value, assume constant cbeh for all planes
    %cellbeh is zeros and ones to determine adhesion or lack thereof
    cell_diam_factor = 1;
    if( (nargin>6) & (length(celldia)>1))
        if (length(celldia)<inp)
           cell_diam_factor=celldia(length(celldia)); 
        else
           cell_diam_factor=celldia(inp);
        end
    end
    
    %now determine adhesion profile
    %if beh>0, use profile from existing distribution, else equal above
    %cell diameter
    
    
    %==============================================
    %PROBABILITY DENSITY FOR ADDITION AND MOVEMENT
    %==============================================
    % the probability density for cell addition is a linear superposition of
    % 1.)    the adhesion profile probDen_AD (high probability in a small range of 
    % permissible distances) and 
    % 2)     the indifference profile probDen_IN (same probability at 
    % all distances, but still exclusion at small distances under the cell 
    % diameter)
    % The relative contributions of the two components are determined by the
    % adhesion behavior parameter, beh. If adhesion is strong (beh=1),
    % addition is only possible at distances like in the first image; if
    % adhesion is weak, these distances are no longer favored as much, or
    % not at all if beh=0. 
    
    % the adhesion profile is stretched, if necessary, by the factor
    % cell_diam_factor which accounts for cell growth
    % if the cells grow significantly, i.e. if the new maximum shifts so
    % far right that it gets to close to the cut-off length, the length has
    % to be extended
    probDen_AD=pdfilter/max(pdfilter);
    current_avcellsize = csiz_max;
    
    if (cell_diam_factor~=1)
        if (cell_diam_factor>1)
            newlength = round(cell_diam_factor * length(pdfilter));
            probDen_AD=zeros(1,newlength);
        end

        llvec = 1:length(probDen_AD);
        probDen_AD(llvec) = pdfilter(round(llvec/cell_diam_factor))/max(pdfilter);
        maxpos_AD = min(find(probDen_AD==max(probDen_AD)));
        current_avcellsize = maxpos_AD;
    end
    
    %indifference probability profile: all distances have the same value,
    %except for the very small distances, which have the size exclusion
    probDen_IN = probDen_AD;
    probDen_IN(current_avcellsize:length(probDen_IN)) = max(probDen_AD);
    probDen_IN = probDen_IN/max(probDen_IN);
         
    %total probability density: linear superposition
    probDen = beh*probDen_AD + (1-beh)*probDen_IN;
    probDen = probDen/max(probDen);
    
    
    
    %====================================================================
    % move points according to allowed probability density
    
    %cycle over all points
    for k=1:nx
        %in this version of the function, cell division has already been 
        %taken care of with the Monte Carlo point addition (see below), so 
        %that the cells only move now
        
        %for every point the length of the movement vector is determined by
        %the average speed, subsequently the angle of movement is tested
        %for compatibility with neighbors
        
        %probDen is adhesion profile for this plane
        %if the cell is an isolated cell, i.e. if the original nearest
        %neighbor distance is nd > 2*csiz+spid, then this restriction is
        %irrelevant
        
        %initialize values
        pointpd = probDen;
        acceptpoint = 0;
        cellIsSingle = 0;
        
        % or cell dist is distance (vector) of the current cell from all other cells
        [orCellDist]=distanceMatrix(mtemp(k,:),mtemp);
        % orCellNNDist is minimum of the vector, i.e. nearest neighbor
        % distance
        orCellNNDist=min(nonzeros(orCellDist));
        
        %orCellNNDist=origNearNeigh(k);
        %if original nearest neighbor distance is larger than ~ twice the 
        %current cellsize, then the cell is 
        %considered single; then its movement is considered to be
        %independent of neighbors, and any movement is therefore accepted
        %
        %if we're in adhesion mode but the cell is too far separate from 
        %other cells, probability is also modified
        if (orCellNNDist > min((2*current_avcellsize+spid),100))
            cellIsSingle = 1;
        end
        
        
        %mpmRest is all other points except the current point
        mpmRestT = mtemp;
        mpmRestT(k,:) = 0;
        mpmRest = [ nonzeros(mpmRestT(:,1)) nonzeros(mpmRestT(:,2)) ];
        
        acceptpoint=0;
        %while the generated point does not yet fulfill the requirements
        %for acceptance, perform loop
        counter = 0;
        while (acceptpoint==0)
        
            %generate vector length
            lenDis = spid+randn(1);
        
            %generate angle
            angle = lastAng(k)+(pi/2)*randn(1);
        
            %resulting endpoint for this movement
            endx = mtemp(k,1) + lenDis * sin(angle);
            endy = mtemp(k,2) + lenDis * cos(angle);
            
            %border correction to yield in-image coordinates
            if endx > matsiz(1), endx = endx-matsiz(1); end
            if endy > matsiz(2), endy = endy-matsiz(2); end
            if endx <= 0, endx = endx+matsiz(1); end
            if endy <= 0, endy = endy+matsiz(2); end
                     
            %calculate nearest neighbor distance to existing cells
            [currCellDist]=distanceMatrix([endx endy],mpmRest);
            currNNdist=min(nonzeros(currCellDist));
            %calculate second-nearest neighbor distance
            currCellDist2 = currCellDist;
            currCellDist2(find(currCellDist2<=currNNdist)) = 0;
            currSNNdist=min(nonzeros(currCellDist2));
            
            
            %determine probability for this distance (if it's in the right
            %ballpark, i.e. no more than 100 pix)
            if (cellIsSingle == 1)
                acceptpoint=1;
            else
                if( currNNdist>min((2*csiz_max+spid),100) )
                    acceptpoint=1;
                else
                    %if second nearest neighbor is too far away to be
                    %relevant
                    if( currSNNdist>min((1.8*currNNdist),100) )
                        %occurence probability for this particular distance
                        currProb = pointpd(ceil(currNNdist));
                    else
                        %occurence probability for this particular distance
                        currProb = pointpd(ceil(currNNdist))*pointpd(ceil(currSNNdist));
                    end
                    %generate random variable between 0 and 1
                    randomVar = rand(1);
                    %do MC test
                    if ( randomVar <  currProb )
                        acceptpoint=1;
                    end
                end% of if
            end % of if
            counter = counter + 1;
            if(counter > 1000)
                acceptpoint=1;
                endx = mtemp(k,1); 
                endy = mtemp(k,2);
            end

        end % of while
        
        %last point is accepted now
        
        mtemp(k,:) = [endx endy];
        lastAng(k) = angle;
        
    end  % of for - cycle over all points in this frame
    
    
    %after moving, decide how many added cells in this frame
    %number of cells is determined by Poisson probability
    expval = gr*nx;
    numaddcell = poissrnd(expval);
    nlen_prev = length(mtemp);
    if (numaddcell>0)
        [mpmMonteCarloAdd]=mpmMCaddsim(mtemp,[matsiz(1) matsiz(2)], numaddcell);
        mtemp = mpmMonteCarloAdd;
        nlen_new=length(mtemp);
        movvec((nlen_prev+1):nlen_new,:)=0;
        lastAng((nlen_prev+1):nlen_new)=0;
    end

    %comment/uncomment the next paragraph to display results of cells in
    %this plane
    plot(mtemp(:,1),mtemp(:,2),'r.');
    axis([-5 matsiz(1)+5 -5 matsiz(2)+5]);
    pause(0.05);
    
    %add mtemp to mpm; account for the fact that the number of points may
    %have changed
    [xm2,ym2]=size(mpm);
    [xmt,ymt]=size(mtemp);
    difx=xmt-xm2;
    if(difx>0)
        addv=zeros(difx,ym2);
        mpm=[mpm; addv];
    end
    mpm=[mpm mtemp]; 
    
end  % of for np number of planes
close(h);

end %of function



function[m2]=distanceMatrix(c1,c2)
%this subfunction makes a neighbour-distance matrix for input matrix m1
%input: c1 (n1 x 2 points) and c2 (n2 x 2 points) matrices containing 
%the x,y coordinates of n1 or n2 points
%output: m2 (n1 x n2) matrix containing the distances of each point in c1 
%from each point in c2
[ncx1,ncy1]=size(c1);
[ncx2,ncy2]=size(c2);
m2=zeros(ncx1,ncx2);
for k=1:ncx1
    for n=1:ncx2
        d=sqrt((c1(k,1)-c2(n,1))^2+(c1(k,2)-c2(n,2))^2);
        m2(k,n)=d;
    end
end
end

