function D=createSparseDistanceMatrix(M,N,threshold,eps)
% createSparseDistanceMatrix calculates and stores all distances <= threshold into a sparse matrix
%
% SYNOPSIS   D=createSparseDistanceMatrix(M,N,threshold,eps)
%
% INPUT      M and N are the matrices containing the set of point coordinates.
%            M and N can represent point positions in 1, 2 and 3D, as follows.
%            
%            In 1D: M=[ y1        and   N=[ y1
%                       y2                  y2
%                       ...                 ... 
%                       ym ]                yn ]
%
%                   Distances: dij = yj-yi
%
%            In 2D:
%                   M=[ y1 x1     and   N=[ y1 x1
%                       y2 x2               y2 x2
%                        ...                 ...
%                       ym xm ]             yn xn ]
%
%                   Distances: dij = sqrt( (yj-yi)^2 + (xj-xi)^2 )
%
%            In 3D:
%                   M=[ y1 x1 z1  and   N=[ y1 x1 z1
%                       y2 x2 z2            y2 x2 z2
%                         ...                ...
%                       ym xm zm ]          yn xn zn ]
%
%                   Distances: dij = sqrt( (yj-yi)^2 + (xj-xi)^2 + (zj-zi)^2 )
%
%            threshold : only the distances dij between the two set of points M and N
%                        which are <= threshold are stored in the (sparse) distance matrix
%
%            eps       : (optional, default value eps = 1e-10)
%                        By definition, sparse matrices contain only non-zero elements.
%                        The sparse matrix returned by createSparseDistanceMatrix contains 
%                        only the distances dij <= threshold. For this reason,
%                        any distance dij > threshold, which is therefore NOT stored in D, is
%                        considered to be equal to 0 by Matlab. 
%                        To be able to distinguish the two cases: dij = 0 from
%                        dij > threshold, all zero distances are replaced in D by a
%                        small number, eps (default eps=1e-10).
%                        Thus the command:
% 
%                           find(D==eps) returns the indices of all zero distances,
%                                        whereas the command:
%                           find(D==0)   returns the indices of all distances > threshold.
% 
% OUTPUT      D        : sparse distance matrix. 
% 
% REMARK      For 1D, both positive and negative distances are returned.
%
% C-MEX file - Aaron Ponti 02/17/2003
