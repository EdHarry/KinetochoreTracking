MATLAB="/Applications/MATLAB_R2012b.app"
Arch=maci64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/Users/edwardharry/.matlab/R2012b"
OPTSFILE_NAME="./mexopts.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for GaussListND_mexCode" > GaussListND_mexCode_mex.mki
echo "CC=$CC" >> GaussListND_mexCode_mex.mki
echo "CFLAGS=$CFLAGS" >> GaussListND_mexCode_mex.mki
echo "CLIBS=$CLIBS" >> GaussListND_mexCode_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> GaussListND_mexCode_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> GaussListND_mexCode_mex.mki
echo "CXX=$CXX" >> GaussListND_mexCode_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> GaussListND_mexCode_mex.mki
echo "CXXLIBS=$CXXLIBS" >> GaussListND_mexCode_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> GaussListND_mexCode_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> GaussListND_mexCode_mex.mki
echo "LD=$LD" >> GaussListND_mexCode_mex.mki
echo "LDFLAGS=$LDFLAGS" >> GaussListND_mexCode_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> GaussListND_mexCode_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> GaussListND_mexCode_mex.mki
echo "Arch=$Arch" >> GaussListND_mexCode_mex.mki
echo OMPFLAGS= >> GaussListND_mexCode_mex.mki
echo OMPLINKFLAGS= >> GaussListND_mexCode_mex.mki
echo "EMC_COMPILER=" >> GaussListND_mexCode_mex.mki
echo "EMC_CONFIG=optim" >> GaussListND_mexCode_mex.mki
"/Applications/MATLAB_R2012b.app/bin/maci64/gmake" -B -f GaussListND_mexCode_mex.mk
