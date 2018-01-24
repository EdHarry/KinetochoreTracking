function [survFunc,pvec] = compSurvivalFunction_general(data1, field1, data2, field2, numBS)
% compSurvivalFunction_general compares the survival functions for
% different conditions
%
% SYNOPSIS [survFunc,pvec] = compSurvivalFunction_general(data1, field1,
% data2, field2)
%
% INPUT     data1    = experiment structure with first data set
%           field1  = field in data1 from which the lifetime histogram is
%                   read, e.g. 'survivalFunction_InRegion'
%           data2   = experiment structure with second data set
%           field2  = field in data1 from which the lifetime histogram is
%                   read, E.G. 'survivalFunction_OutRegion'
%           numBS  = number of bootstrap runs (OPTIONAL)
%                    default = 2000
%                   
%
% OUTPUT    survFunc = survival function 
%                       (first row DATA1, second row Data2)
%           pvec    = vector of p-values;
%                     1: KS-test on difference between data1 and data2
%                     2: 95% threshold for KS-test on bootstrapped
%                     subsamples of data1 vs bootstrapped subsamples of
%                     data2
%                     3. fraction of non-significant boostrap KS-test
%                     p-values 
%
% NOTE: The current version of the function assumes that ALL MOVIES in the
% data structure are acquired at the same framerate, or that the user
% chooses to treat them as being the same, and that the survival functions
% are normalized to the value for the lifetime 1 frame.
%
%
% last modified DATE: 28-Aug-2008 (Dinah)


% averaging can only be performed up until the minimum common length of all 
% movies, so we determine the shortest movie length in the structure
n1 = length(data1);
n2 = length(data2);

for i=1:n1
    mlvec1(i) = data1(i).movieLength;
end
minlen1 = min(mlvec1);

for i=1:n2
    mlvec2(i) = data2(i).movieLength;
end
minlen2 = min(mlvec2);

minlen = min(minlen1,minlen2);


histmat1 = nan*zeros(length(data1),minlen);
histmat2 = nan*zeros(length(data2), minlen);

% read data for each movie 
for k=1:n1
    
    if isfield(data1,field1)
        survFun = getfield(data1(k), field1);
        histmat1(k,:) = survFun(1:minlen);
    else
        error(['function requires a structure field called ',field1]);
    end          
end

% read data for each movie 
for k=1:n2
    
    if isfield(data2,field2)
        survFun = getfield(data2(k), field2);
        histmat2(k,:) = survFun(1:minlen);
    else
        error(['function requires a structure field called ',field2]);
    end          
end


% average survival functions
survFun1_SUM = nansum(histmat1,1);
survFun2_SUM = nansum(histmat2,1);
survFun1_AVE = survFun1_SUM/n1;
survFun2_AVE = survFun2_SUM/n2;

% define output
survFunc(1,:) = survFun1_AVE;
survFunc(2,:) = survFun2_AVE;


% convert to distributions
df1 = convSurvivalFunction2dist(survFun1_SUM);
df2 = convSurvivalFunction2dist(survFun2_SUM);
nd1 = length(df1);
nd2 = length(df2);

% compare averages
[H,pval_av] = kstest2(df1,df2);




%=================================
%% bootstrap

% number of bootstrap runs
if nargin>4
    nbs = numBS;
else
    nbs = 2000;
end


% bootstrap first data set
for b=1:nbs
    
    fprintf('bootstrap %03d',round(100*(b/nbs)));
    
    %=================================
    % bootstrap from reshuffled movies
    
    % position vector
    pos1 = randsample(n1,n1,true); 
    % local bootstrap set from normalized matrix
    cmat1 = histmat1(pos1,:);  
    bootstrapSUM1 = nansum(cmat1,1); 
    df1_sub = convSurvivalFunction2dist(bootstrapSUM1);

    % position vector
    if n2==n1
        pos2 = pos1;
    else
        pos2 = randsample(n2,n2,true); 
    end
    % local bootstrap set from normalized matrix
    cmat2 = histmat2(pos2,:); 
    bootstrapSUM2 = nansum(cmat2,1); 
    df2_sub = convSurvivalFunction2dist(bootstrapSUM2);
    
    %===========================================
    % comp between subsample 1 and total sample 2
    %[H,pval_bs1] = kstest2(df1_sub,df2);
    [H,pval_bs1] = kstest2(df1_sub,df2_sub);
    bootstrap1_pval(b) = pval_bs1;
    
     
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b');
    
end

fprintf('\n');


% sorted bootstrap confidence intervals
sort_bs1 = sort(bootstrap1_pval);

% 95% of all bootstrap comparison p-values are smaller than...
cflevel = sort_bs1(round(0.95*nbs));

% the fraction of bootstrap comparison p-values that are non-significant
% (i.e. larger than 0.05) is...
plevel = length(find(sort_bs1>0.05))/nbs;



% output:
pvec(1) = pval_av; % KS-test p-value of all data1 vs all data 2
pvec(2) = cflevel; % KS-test p-value below which 95% of bootstrap tests are located
pvec(3) = plevel;  % percentage of non-significant boostrap test results




%==========================================================================
%% display results

figure; hold on;
cdfplot(df1);
cdfplot(df2);
h = findobj(gca,'type','line');
set(h(1),'linestyle',':','color','r')

xlabel('lifetime (frames)');
ylabel('cumulative fraction');
legend('data1','data2');
axis([0 minlen -0.01 1.01]);
box on


textstr1 = ['ave1 vs ave2: KS-test p=',num2str(pval_av)];

textstr2 = ['bootstrap1 vs bootstrap2 95% p<',num2str(cflevel)];

if length(find(sort_bs1>0.05)) == 0
    textstr3 = ['fraction of non-sign. bootstraps < ',num2str(1/nbs)];
else
    textstr3 = ['fraction of non-sign. bootstraps = ',num2str(plevel)];
end

text(60,0.5,textstr1);
text(60,0.4,textstr2);
text(60,0.3,textstr3);

grid off



end % of function




%=========================================================================
%
%                            SUBFUNCTION



function [df] = convSurvivalFunction2dist(sf)
% convSurvivalFunction2dist converts a survival function back into a
% distribution function
%
% SYNOPSIS [df] = convSurvivalFunction2dist(sf)
%
% INPUT     sf    = survival function, presumed at equal spacing 
%
% OUTPUT    df    = distribution function
%
% NOTE: Advantage of this function is that survival functions can be summed
% and cropped until appropriate time points, that that the reconversion
% gets rid of the problem of unequal movie lengths


% convert to histogram
hf = abs(diff(sf));

df = zeros(sum(hf),1);
ct = 0;

for i=1:length(hf)
    numentry = round(hf(i));
    %cv = i + zeros(numentry,1);
    df(ct+1:ct+numentry) = i; 
    ct = ct+numentry;
end % of for

end % of subfunction




    
