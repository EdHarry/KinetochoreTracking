function interField=vecFldInterpAnisoB(M, tgtPoints, sigmaU, sigmaV, maxIteration, polygon)
% vecFldInterpAnisoB computes the vector field on a user-specified point set based on vectors
% on a given point set. Initial direction of the flow vector at each given
% point is computed using an isotropic filter. Then this initial vector is
% used for anisotropic filtering at this point. This process is iterated
% until either convergence (in 2-norm) is reached or the maximum number
% iteration is reached. This function do not require users to provide
% initial external flow field. However, if this is desired, function
% "vecFldInterpAnisoA" can be used. 
%
% Interpolation is based on an anisotropic Gaussain probability density function of the following form:
%  
% p(u, v) =  1 / (2 * pi * sigmaU * sigmaV) * exp(-0.5 * (u / sigmaU)^2 -
%            0.5 * (v / sigmaV)^2);
%
% SYNOPSIS   interpFiled = vecFldInterpAnisoB(M, tgtPoints, sigmaU, sigmaV, maxIteration, polygon)
%
% INPUT              M  :     known vector field, stored in a (nx4)-matrix of the form [y0 x0 y x]n
%                             (where (y0,x0) is the base and (y,x) is the tip of
%                             the vector).
%
%            tgtPoints  :     points at which vector field needs to be calculated.
%                             It is stored in a (mx2)-matrix of the form [yg xg]m.In particular, 
%                             this point set can be the same as the given point set in M.
%
%       sigmaU, sigmaV  :     Gauss distribution parameters for calculating
%                             the probability function. Currently, it is
%                             assumed that they are contant over the entire
%                             field. sigmaU is defined along the flow
%                             direction. sigmaV is defined perpendicular to
%                             the flow direction.
%                             NOTICE: these two parameters will also be used to determine a search region, based 
%                             essentially on the three-sigma rule. 
%
%         maxIteration  :     maximum number of iterations. This number may not be reached as the calculation
%                             of flow vector often converges faster. 
% 
%              polygon  :     (optional - pass polygon=[] to disable). The interpolated vector
%                             can be cropped to remove vectors outside a given region of interest.
%                             To create the polygon use the functions ROIPOLY or
%                             GETLINE. These functions return the polygon vertices
%                             stored in two vectors y and x. Set polygon=[y x] to
%                             use with vectorFieldInterp.
%
% OUTPUT    interField  :     interpolated vector field for the points given in tgtPoints.

% DEPENDENCES                 This function calls vecFldInterpAnisoA
%
% REMARKS                     This implementation is based on the function and
%                             interface of "vectorFieldInterp" created by Aaron Ponti,
%                             11/18/2002. 
%
% AUTHOR                      Ge Yang 
%                             geyang@scripps.edu
%                             Laboratory for Computational Cell Biology
%                             The Scripps Research Institute
%
% DATE OF CREATION            March 10, 2004
% 
% CHANGES
%                             ver 0.3: The code is reorganized for better
%                             readability and initial release.
%
%
%                             ver 0.4: use MEX function to compute distance                                 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
% Vector base positions
srcPoints = M(:, 1:2); % These will be the centers of the probability ellipsoids.
V=[M(:, 3) - M(:, 1) M(:, 4) - M(:, 2)];
SCALEFACTOR = 3;   % use 3-sigma rule

[srcLength temp0] = size(srcPoints);
[tgtLength temp0] = size(tgtPoints);

r = 9 * (sigmaU^2 + sigmaV^2); % choose a relatively large one for strong smoothing effect
sigma = sigmaU^2 + sigmaV^2;     % still need a Gauss kernel, though isotropic in this case.

% Step 1: Generate the initial flow field using an isotropic filter. 

% First, call the MEX function to compute distance. This will save
% computation time.

D = createDistanceMatrix(srcPoints, tgtPoints);

prevFlowVecList = zeros(tgtLength, 2);
for i = 1 : tgtLength
    weight_sum = 0;
    vecsum = [0 0];
    for j = 1 : srcLength
        temp = D(j, i);  % the Euclidean distance
        if temp <= r  % actually comparing distance square
            weight = exp(-0.5 * temp / sigma); % will normalize later. So no need to compute leading coefficient
            vecsum = vecsum + weight * V(j, :);
            weight_sum = weight_sum + weight;
        end
    end
    prevFlowVecList(i, :) = vecsum / weight_sum;   % This is our initial vector field
end
clear D;
% Step 2: Call vecFldInterpAnisoA iterative to compute the vector field

error = 10^6; % just give a large number. value not meaningful
threshold = 1e-2; 

i = 0;
while (1)
    newFlowVecList = vecFldInterpAnisoA(M, tgtPoints, prevFlowVecList, sigmaU, sigmaV);
    error = sqrt(sum(newFlowVecList - prevFlowVecList).^2) / tgtLength;
    if (error < threshold)
        interField = newFlowVecList;
        fprintf('Computation converged in flow field calculation. Now return to the calling point. \n');
        return;
    else
        i = i + 1;
        if (i >= maxIteration)
            fprintf('WARNING: Maximum number reached in iterative calculation of vector field. Now return to the calling point. \n');
            interField = newFlowVecList;
            return;
        end
    end
    fprintf('Flow field calculation: iteration # %d\n', i);
    prevFlowVecList = newFlowVecList;
end



