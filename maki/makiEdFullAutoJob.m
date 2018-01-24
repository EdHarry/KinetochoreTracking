function makiEdFullAutoJob(jobVector)
% written by Ed 14/04/10, to automattically setup a job, analyse movies and
% produce output (uses all default settings for sister fitting etc and asumes movies are already corectly cropped)

tic

[fileName,dir2SaveRes] = uiputfile('makiAnalysis_*');

[job,basePath]=makiMakeJob_ed('MCAINSH',jobVector,[],0);







makiMovieAnalysis('MCAINSH',job)

analysisStruct = makiCollectMovies_ed('MCAINSH',basePath,fileName,dir2SaveRes);

try
analysisStruct = makiSepDispSpaceTimeAnalysis('MCAINSH',analysisStruct,1);
end

try
analysisStruct = makiSisterMotionCoupling('MCAINSH',analysisStruct,1);
end

try
analysisStruct = makiSisterConnectionAnalysis('MCAINSH',analysisStruct,1);
end

try
analysisStruct = makiDiffusionAnalysis ('MCAINSH',analysisStruct,1);
end

toc

end