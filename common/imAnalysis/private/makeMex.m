%% compile commands for all mex files in this directory
clear mex;
mex mexGetLine.c tbx.lib -f mexOpts.bat
%mex mexNoiseEstim.c ../../mexCalls/mexDataConverters.c tbx.lib -f mexOpts.bat
mex -g mexNoiseEstim.c mexUtils.lib tbx.lib -f mexOpts.bat 
mex -g mexRayleighMode.c mexUtils.lib tbx.lib -f mexOpts.bat 
%mex mexLineDetect.c ../../mexCalls/mexDataConverters.c tbx.lib -f mexOpts.bat
mex -g mexLineDetect.c mexUtils.lib tbx.lib -f mexOpts.bat