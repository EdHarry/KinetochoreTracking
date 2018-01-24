MATLAB="/Applications/MATLAB_R2013b.app"
Arch=maci64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/Users/edwardharry/.matlab/R2013b"
OPTSFILE_NAME="./mexopts.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for mmfMex" > mmfMex_mex.mki
echo "CC=$CC" >> mmfMex_mex.mki
echo "CFLAGS=$CFLAGS" >> mmfMex_mex.mki
echo "CLIBS=$CLIBS" >> mmfMex_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> mmfMex_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> mmfMex_mex.mki
echo "CXX=$CXX" >> mmfMex_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> mmfMex_mex.mki
echo "CXXLIBS=$CXXLIBS" >> mmfMex_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> mmfMex_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> mmfMex_mex.mki
echo "LD=$LD" >> mmfMex_mex.mki
echo "LDFLAGS=$LDFLAGS" >> mmfMex_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> mmfMex_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> mmfMex_mex.mki
echo "Arch=$Arch" >> mmfMex_mex.mki
echo OMPFLAGS= >> mmfMex_mex.mki
echo OMPLINKFLAGS= >> mmfMex_mex.mki
echo "EMC_COMPILER=" >> mmfMex_mex.mki
echo "EMC_CONFIG=optim" >> mmfMex_mex.mki
"/Applications/MATLAB_R2013b.app/bin/maci64/gmake" -B -f mmfMex_mex.mk
