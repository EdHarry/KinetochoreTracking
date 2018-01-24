function aObjectId = StartImarisInstance_v73
% EHarry Nov 2011
unix('/Users/Ed/Documents/MATLAB/maki/imarisICE/startImaris_v73.sh');
!/Applications/Imaris\ 7.3.0.app/Contents/MacOS/ImarisServerIce &
aObjectId = 0;
end