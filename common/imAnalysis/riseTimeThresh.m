function[t, lev]=riseTimeThresh(fptrace,norm);
% function[t, lev]=triseFirstpoint(fptrace,norm);
% calculates thresholded risetime of continuously rising time series
% 
% INPUT: fptrace  = continuous trace time series (e.g. velocity or scattering
%                   parameter) - this trace has to be normalized to values
%                   between 1 and 2 for this function to yield useful
%                   values 
%        norm = normalization time point (frame in which drug is added); by
%               definition, changes on the left of this time point are not
%               considered to be persistent rises
%        
%       IMPORTANT: consider the possibility that due to windowing in e.g.
%        the scatter parameter determination, the time point of drug 
%        addition may actually have been reduced in rleation to the original
% 
% OUTPUT: t = rise time vector with number of frames per level; 
%               although the rise time level trace considers 10 levels as 
%               specified below, the time vector only contains the rise times 
%               for the first 6 level steps - of course this can be changed
%               according to the user's needs
%        lev  = copy of fptrace with discrete levels
%
%last changed June 09, 2005 by Dinah Loerke


%==============================================================
%
%   define levels for rise time
%
%   in the current implementation, the function assumes an input
%   function that is normalized to 1, and rises to ~ 2 
%
%=============================================================
levelt1=1.1;
levelt2=1.2;
levelt3=1.3;
levelt4=1.4;
levelt5=1.5;
levelt6=1.6;
levelt7=1.7;
levelt8=1.8;
levelt9=1.9;
levelt10=2.0;

len=length(fptrace);

%==============================
% nan values of the original trace are filled using adjacent points
%===============================

for i=1:len
    if(isnan(fptrace(i))==1)
        a=min(i,abs(i-5));
        b=min(i+5,len);
        fptrace(i)=nanmean(fptrace(a:b));
    end
end
wtrace=fptrace;

%==================================
%   original trace (fptrace) is filtered, yielding wtrace
%==================================

windowsize = 21;
orfiltertrace=filter([1 1 1 ]/3,1,fptrace);
filtertrace=orfiltertrace((2):len);
wtrace=filtertrace;

%==================================
%   thresholded copy is calculated (leveltrace)
%==================================

leveltrace=wtrace;
leveltrace(wtrace<=levelt1)=1;
leveltrace(wtrace>levelt1)=levelt1;
leveltrace(wtrace>levelt2)=levelt2;
leveltrace(wtrace>levelt3)=levelt3;
leveltrace(wtrace>levelt4)=levelt4;
leveltrace(wtrace>levelt5)=levelt5;
leveltrace(wtrace>levelt6)=levelt6;
leveltrace(wtrace>levelt7)=levelt7;
leveltrace(wtrace>levelt8)=levelt8;
leveltrace(wtrace>levelt9)=levelt9;
leveltrace(wtrace>levelt10)=levelt10;
leveltrace(1:norm)=1;


%======================================
%   for each level change in the above calculated leveltrace, look if the
%   jump is significant, i.e. persistent over time
%========================================

dlevel = diff(leveltrace);
dlevelc = dlevel;
leveltracec = leveltrace;
%loop over all time points
for i=1:(length(dlevel)-windowsize)
    %check if there is a level increase or decrease at this point (i.e. a
    %change not already previously erased, so use dlevelc instead of dlevel)
    if(dlevelc(i)~=0)
        %check if this change is at some point later compensated, so that the
        %function returns to its original level
        %start checking from smaller to larger windowsizes; if smaller
        %windowsize is applicable, do not go checking for larger sizes
        %(meaning e.g. that if there are several transient maxima or minima
        %contained in windowsize, choose only the first one, and consider
        %the later one separately, with their own winsowsize "future"
        for k=1:windowsize
            %however, due to the typical shape of the function, we consider 
            %an increase significant if there's a double increase in the
            %window; in that case, we exit the loop
            if( sum(dlevel(i:i+k))>= 1.5*abs(dlevelc(i)))
                break                
            end
            %if there's not double increase but a compenstation, we 
            %consider this maximum or minimum not significant and erase it 
            if( sum(dlevel(i:i+k))==0 )
                setfwinsiz = k;
                dlevelc(i:i+setfwinsiz)=0;
                break
            end %of if
        end % of for
     end %of if levelchange occurs
end % of for

%========================
%   modify the calculated leveltrace using the results of the persistence
%   test above
%=======================

for i=2:length(dlevel)
    leveltracec(i) = leveltracec(i-1)+dlevelc(i-1);
end

%=======================================================
%    comment this paragraph is you want to forego the plotting
%=======================================================

plot(fptrace,'b.')
hold on
axis([0 190 0.9 2.05]);
plot(filtertrace,'go')
plot(leveltrace,'b-')
plot(leveltracec,'r-')
hold off

%leveltracec is shorter than original because of filtering (-2)
evaltrace = leveltracec (norm:(len-2));
indexw=1:(len-norm-1);

%==========================
%    now determine risetimes in modified leveltrace
%
%   IMPORTANT: if you want to include additional rise time levels, you have
%   to add the appropriate lines here for their calculation, e.g. t7=...,
%   and then add t7 to the vector in line 181 : t=[t1 t2 t3 t4 t5 t6 t7 t8 ...];
%===========================
t1=min(indexw(evaltrace>=levelt1));
t2=min(indexw(evaltrace>=levelt2))-(t1);
t3=min(indexw(evaltrace>=levelt3))-(t1+t2);
t4=min(indexw(evaltrace>=levelt4))-(t1+t2+t3);
t5=min(indexw(evaltrace>=levelt5))-(t1+t2+t3+t4);
t6=min(indexw(evaltrace>=levelt6))-(t1+t2+t3+t4+t5);
if(isempty(t1))
    t1=len;
end
if(isempty(t2))
    t2=0;
end
if(isempty(t3))
    t3=0;
end
if(isempty(t4))
    t4=0;
end
if(isempty(t5))
    t5=0;
end
if(isempty(t6))
    t6=0;
end

%==================
%    return values
%==================
t=[t1 t2 t3 t4 t5 t6];
lev = [ones(norm,1)' evaltrace];
