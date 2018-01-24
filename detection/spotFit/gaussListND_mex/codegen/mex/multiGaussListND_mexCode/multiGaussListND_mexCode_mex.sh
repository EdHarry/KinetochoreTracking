MATLAB="/Applications/MATLAB_R2012b.app"
Arch=maci64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/Users/edwardharry/.matlab/R2012b"
OPTSFILE_NAME="./mexopts.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for multiGaussListND_mexCode" > multiGaussListND_mexCode_mex.mki
echo "CC=$CC" >> multiGaussListND_mexCode_mex.mki
echo "CFLAGS=$CFLAGS" >> multiGaussListND_mexCode_mex.mki
echo "CLIBS=$CLIBS" >> multiGaussListND_mexCode_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> multiGaussListND_mexCode_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> multiGaussListND_mexCode_mex.mki
echo "CXX=$CXX" >> multiGaussListND_mexCode_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> multiGaussListND_mexCode_mex.mki
echo "CXXLIBS=$CXXLIBS" >> multiGaussListND_mexCode_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> multiGaussListND_mexCode_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> multiGaussListND_mexCode_mex.mki
echo "LD=$LD" >> multiGaussListND_mexCode_mex.mki
echo "LDFLAGS=$LDFLAGS" >> multiGaussListND_mexCode_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> multiGaussListND_mexCode_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> multiGaussListND_mexCode_mex.mki
echo "Arch=$Arch" >> multiGaussListND_mexCode_mex.mki
echo OMPFLAGS= >> multiGaussListND_mexCode_mex.mki
echo OMPLINKFLAGS= >> multiGaussListND_mexCode_mex.mki
echo "EMC_COMPILER=" >> multiGaussListND_mexCode_mex.mki
echo "EMC_CONFIG=optim" >> multiGaussListND_mexCode_mex.mki
"/Applications/MATLAB_R2012b.app/bin/maci64/gmake" -B -f multiGaussListND_mexCode_mex.mk
