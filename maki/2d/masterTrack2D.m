function masterTrack2D

% EHarry Oct 2011

%% ORIGIANL HEADER
% % % written by Ed 14/04/10, to automattically run a job on parallel cores
% % % produce output (uses all default settings for sister fitting etc and asumes movies are already corectly cropped)

jobVector=[1 2 3 5 6 7 8 9];

[fileName,dir2SaveRes] = uiputfile('makiAnalysis_*'); %select location to save results

[job,basePath]=makiMakeJob_ed('MCAINSH',jobVector,[],0); %makimakejobed returns the basepath aswell as the job struct, it also automatically uses default tracking parameters, change last 0 to 1 to make fuction ask for parameters

tic %start timer

no_jobs=length(job);

%no_workers = no_jobs;


batch_jobs=struct('jobs',0);



warning off
for i=1:no_jobs
    
    s= ['makiMovieAnalysis2D(job(' int2str(i) '),1);']; %makimovie analysis_ed runs default setings automatically, and only needs the job struct passed in
    batch_jobs(i).jobs=batch(s); %the batch command sends jobs to the matlab queue
    
    
end

warning on





for i=1:no_jobs
    x=batch_jobs(i).jobs; %waits on all the jobs to finish
    wait(x);
end

try
    analysisStruct = makiCollectMovies_ed('MCAINSH',basePath,fileName,dir2SaveRes); %makicollectmovies_ed automatically collects movies in the basepath
    
    try
        analysisStruct = makiSepDispSpaceTimeAnalysis('MCAINSH',analysisStruct,1); %tries to run all analysis on the movies
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
end

toc %end timer

end