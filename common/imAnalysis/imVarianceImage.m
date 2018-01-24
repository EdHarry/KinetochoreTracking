function [sigmaSquare] = imVarianceImage (varargin)
% VarianceImage takes an input image and calculates the variance of every pixel.
%
% SYNOPSIS       [outImage] = imVarianceImage (varargin)
%
% INPUT          inputImage: the original greylevel image including background
%                W : odd integer denoting the width of the neighbourhood
%
% OUTPUT         sigmaSquare: a normalized image where every pixel is represented by the variance
%                             of the intensity of the corresponding input image pixel
%
% DEPENDENCIES   imVarianceImage uses { nothing }
%                                  
%                imVarianceImage is used by { ptGetProcessedImage }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Apr 04          Initial release
% Andre Kerstens        May 04          Normalized the output image  

% Do some error checking on the input
if nargin < 2
   error('The input image and neighbourhood W (odd int) have to be provided. See help imVarianceImage.');
end

inputImage = varargin{1};
W = varargin{2};

% Calculate InputImage^2
inputSquare = inputImage.^2;

% Generate a mask out of the neighbourhood width
mask = ones (W) / W^2;

% Convolve the input image with the mask
inputImageConv = conv2 (inputImage, mask, 'same');

% Convolve the squared input image using the same mask
inputSquareConv = conv2 (inputSquare, mask, 'same');

% Subtract the convolution of the squared input from the squared convolution of the input
% to get the variance
sigmaSquare = inputSquareConv - inputImageConv.^2;

% As a last step normalize the variance matrix
minSigmaSquare = min (min (sigmaSquare));
maxSigmaSquare = max (max (sigmaSquare));
sigmaSquare = (sigmaSquare - minSigmaSquare) / (maxSigmaSquare - minSigmaSquare);
