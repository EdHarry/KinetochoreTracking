function aObjectId = StartImarisInstance
% EHarry Nov 2011
unix('/Users/Ed/Documents/MATLAB/maki/imarisICE/startImaris.sh');
!/Applications/Imaris\ 7.4.0.app/Contents/MacOS/ImarisServerIce &
aObjectId = 0;
end