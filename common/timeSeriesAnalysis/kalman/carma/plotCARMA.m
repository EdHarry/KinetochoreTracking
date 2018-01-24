function [] = plotCARMA(TRAJ)
%Plots the multivariate time series TRAJ in a single figure in pretty
%colors

%Uses longest dimension as time and next as time series index


figure(1)
%Check that matrix is 2D
if ndims(TRAJ) > 2
    disp('Input time series matrix should be 2D');
    return
end

%Make sure that time is 1st dimension
if size(TRAJ,1) < size(TRAJ,2)
    TRAJ = TRAJ'
end

%Determine number of timeseries
nNode = size(TRAJ,2);
subplot(nNode,1,1)

%Plot the time series
for j = 1:1:nNode
   figure(1);
   subplot(nNode,1,j)
   plot(TRAJ(:,j),'color',[(1-(j/nNode)) 0 (j/nNode) ]);
   xlabel('Time');
   ylabel('Amplitude (A.U.)');   
end

