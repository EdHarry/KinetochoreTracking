function u = nonLinearFit(x,y,u0)
% NONLINEARFIT uses non-linear least squares to fit input data. 
% The model has to be specified in the subfunction 'fittingFcn'. 
% Type 'edit nonLinearFit' and scroll down to enter the fittingFcn.
% To load data from an Excel-sheet, use 'A=xlsread('C:\data...\myFile.xls')'.
% If the data is stored as [x, y], x=A(:,1), y=A(:,2).
%
% SYNOPSIS u = nonLinearFit(x,y,u0)
% 
% INPUT    x : x-values of the data points
%          y : y-values of the data points
%          u0: initial guess for the unknowns. The length of vector u0 has
%               match the number of unknowns in your problem. If you have
%               no idea of the value of the unknowns, take ones.
%
% OUTPUT   u : fitted unknowns
%  
% REMARKS
%  The program will also plot x,y with the fitted curve, and the residuals.
%  Of course, the program will also perform linear least squares fits.
%
% c: 10/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=============
% TEST INPUT
%=============

% check number of input arguments
if nargin < 3 || isempty(x) || isempty(y) || isempty(u0)
    error('nonLinearFit needs three non-empty input arguments');
end

% make sure x and y have the same size
x=x(:);
y=y(:);

if length(x) ~= length(y)
    error('x and y have must have the same number of elements');
end

% make sure that we have enough data points
if length(u0) > length(x)-1
    error('there are not enough data points to fit so many unknowns')
end

%===============


%===========================
% CALL THE FITTING FUNCTION
%===========================

% echo off
options = optimset('Display','off');
u = lsqnonlin(@fittingFcn, u0, [], [], options, x, y);

% plot function
residuals = fittingFcn(u,x,y);
figure('Name','Data and fit');
plot(x,y,'.k',x,y-residuals,'-r');

% plot residuals
figure('Name','Residuals')
stem(x,residuals)

%===========================





%=============================
% FITTING FUNCTION
%=============================
function residuals = fittingFcn(u, x, y)
% write here the residuals, i.e. y-f(x,u), where u are the unknowns that
% you want to fit. If there are many unknowns, u is a vector, so the first
% element will be u(1), the second u(2) etc. Don't forget that x and y are
% vectors, too, so for multiplications or powers you have to use
% element-wise operations (e.g. .*, .^)!
% Lines starting with a '%' sign are regarded as comments, so if you have
% to use another function, but not retype the current function the next
% time you use it, just comment it out.

% ---- calculation of the function values f(x,u)

% exponential fit
fcn = u(1) * exp(x/u(2));


% % linear fit (just for fun)
% fcn = u(1) + u(2) * x + u(3)*x.^2;


% --- calculation of the residuals - this line has to be at the very end!
residuals = y - fcn;