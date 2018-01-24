function master_track_poleOME(job,jobVector,updateJob,askDecon)
% written by EHarry 14/04/10, to automattically run a job on parallel cores
% produce output (uses all default settings for sister fitting etc and asumes movies are already corectly cropped)

if nargin < 1
    job = [];
end

if nargin < 2 || strcmp(jobVector,'all')
    jobVector=[1 2 3 5 6 7 8 9 10 11 12 13 14 15 16 17];
end

if nargin < 3 || isempty(updateJob)
    updateJob = 1;
end

if nargin < 4 || isempty(askDecon)
    askDecon = 1;
end
    
%[fileName,dir2SaveRes] = uiputfile('makiAnalysis_*'); %select location to save results

if updateJob || isempty(job)
    job = makiMakeJob_poleOME('MCAINSH',jobVector,job,0,[],askDecon); %makimakejobed returns the basepath aswell as the job struct, it also automatically uses default tracking parameters, change last 0 to 1 to make fuction ask for parameters
end

%tic %start timer

no_jobs=length(job);

%no_workers = no_jobs;


batch_jobs=struct('jobs',0);



warning off
for i=1:no_jobs

 s= ['makiMovieAnalysis_poleOME_auto(0,job(' int2str(i) '),3);']; %makimovie analysis_ed runs default setings automatically, and only needs the job struct passed in 
  batch_jobs(i).jobs=batch(s); %the batch command sends jobs to the matlab queue


end

warning on



progressText(0,'\nAnalysing Movies\n');

for i=1:no_jobs
    x=batch_jobs(i).jobs; %waits on all the jobs to finish
    wait(x);
    progressText(i/(no_jobs),sprintf('\nAnalysed %i Movies\n',i));
end


% analysisStruct = makiCollectMovies_ed('MCAINSH',basePath,fileName,dir2SaveRes); %makicollectmovies_ed automatically collects movies in the basepath
% 
% try
% analysisStruct = makiSepDispSpaceTimeAnalysis('MCAINSH',analysisStruct,1); %tries to run all analysis on the movies
% end
% 
% try
% analysisStruct = makiSisterMotionCoupling('MCAINSH',analysisStruct,1);
% end
% 
% try
% analysisStruct = makiSisterConnectionAnalysis('MCAINSH',analysisStruct,1);
% end
% 
% try
% analysisStruct = makiDiffusionAnalysis ('MCAINSH',analysisStruct,1);
% end
% 
% 
% toc %end timer

end