MATLAB="/Applications/MATLAB_R2013b.app"
Arch=maci64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/Users/edwardharry/.matlab/R2013b"
OPTSFILE_NAME="./mexopts.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for fitNGaussians3D_mexCode" > fitNGaussians3D_mexCode_mex.mki
echo "CC=$CC" >> fitNGaussians3D_mexCode_mex.mki
echo "CFLAGS=$CFLAGS" >> fitNGaussians3D_mexCode_mex.mki
echo "CLIBS=$CLIBS" >> fitNGaussians3D_mexCode_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> fitNGaussians3D_mexCode_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> fitNGaussians3D_mexCode_mex.mki
echo "CXX=$CXX" >> fitNGaussians3D_mexCode_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> fitNGaussians3D_mexCode_mex.mki
echo "CXXLIBS=$CXXLIBS" >> fitNGaussians3D_mexCode_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> fitNGaussians3D_mexCode_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> fitNGaussians3D_mexCode_mex.mki
echo "LD=$LD" >> fitNGaussians3D_mexCode_mex.mki
echo "LDFLAGS=$LDFLAGS" >> fitNGaussians3D_mexCode_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> fitNGaussians3D_mexCode_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> fitNGaussians3D_mexCode_mex.mki
echo "Arch=$Arch" >> fitNGaussians3D_mexCode_mex.mki
echo OMPFLAGS= >> fitNGaussians3D_mexCode_mex.mki
echo OMPLINKFLAGS= >> fitNGaussians3D_mexCode_mex.mki
echo "EMC_COMPILER=" >> fitNGaussians3D_mexCode_mex.mki
echo "EMC_CONFIG=optim" >> fitNGaussians3D_mexCode_mex.mki
"/Applications/MATLAB_R2013b.app/bin/maci64/gmake" -B -f fitNGaussians3D_mexCode_mex.mk
