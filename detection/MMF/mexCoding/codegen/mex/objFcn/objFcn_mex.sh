MATLAB="/gpfs/hpcwarwick/matlab/r2011a"
Arch=glnxa64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/home/sysbio/msrfby/.matlab/R2011a"
OPTSFILE_NAME="./mexopts.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for objFcn" > objFcn_mex.mki
echo "CC=$CC" >> objFcn_mex.mki
echo "CFLAGS=$CFLAGS" >> objFcn_mex.mki
echo "CLIBS=$CLIBS" >> objFcn_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> objFcn_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> objFcn_mex.mki
echo "CXX=$CXX" >> objFcn_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> objFcn_mex.mki
echo "CXXLIBS=$CXXLIBS" >> objFcn_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> objFcn_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> objFcn_mex.mki
echo "LD=$LD" >> objFcn_mex.mki
echo "LDFLAGS=$LDFLAGS" >> objFcn_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> objFcn_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> objFcn_mex.mki
echo "Arch=$Arch" >> objFcn_mex.mki
echo OMPFLAGS= >> objFcn_mex.mki
echo OMPLINKFLAGS= >> objFcn_mex.mki
echo "EMC_COMPILER=unix" >> objFcn_mex.mki
echo "EMC_CONFIG=optim" >> objFcn_mex.mki
"/gpfs/hpcwarwick/matlab/r2011a/bin/glnxa64/gmake" -B -f objFcn_mex.mk
