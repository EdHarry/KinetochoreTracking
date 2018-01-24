function h= mywaitbar(op,h,maxval,name);
%wrapper for the built-in function waitbar with cancel-button and ETA
%
%SYNOPSIS h = mywaitbar(op,h,maxval,name)
%
%INPUT   op     operation that is carried out
%                   -1 just draw again, no update
%                   -2 abort
%                   otherwise specify how much of the job is done
%                   (e.g. loopCounter/maxValueOfCounter)
%        h      waitbar-handle (pass empty to launch waitbar)
%        maxval not needed anymore
%        name   title of waitbar figure (optional)
%
%OUTPUT  h      waitbar-handle
%
%---------How to use: 
%
%outside of the loop call 
%waitbarHandle=mywaitbar(0,[],[],name)
%
%within the loop with a controlVariable going from 1:maxval
%to update and to be able to benefit most from cancel-button:
%
%try
%begin loop
%do calculations
%mywaitbar(controlVariable/maxValueOfControlVariable,waitbarHandle,[]);
%end loop
%catch
%if findstr(lasterr,['Error using ==> get',char(10),'Invalid handle'])
%--insert here the action you want to happen on pressing cancel
%else
%rethrow(lasterror) %or any other action to be executed on an error within the loop
%end
%
%close(waitbarHandle)
%
%---------------------
%
%c: dT 01/03

if isempty(h)
    h = waitbar(0,'Please wait...');
    set(h,'Units','pixels');
    pb=uicontrol(h,'Style','PushButton','String','Cancel','Position',[160 0 40 20],'Callback','mywaitbar(-2,0,0);');
    if nargin==4
        set(h,'Name',name);
    end        
    tic;
else
    t=toc;
    if t>0.1  %min update time
        switch op
        case -1 %just figure refresh
            drawnow;
        case -2 %break
            h=findobj('Tag','TMWWaitbar');
            close(h)
            error('user aborted process...')
        otherwise  %update pos or just figure refresh
            timeSeries=get(h,'UserData');
            timeSeries=[timeSeries toc]; %get time since last call in loop
            set(h,'UserData',timeSeries);
            %calculate ETA: avg time of one step * expected # of remaining
            %steps (= current# of steps/current percentage*remaining percentage)
            remainingTime = median(timeSeries) * length(timeSeries) * (1-op)/(op+0.00001); %add a little for op==0
            waitbar(op,h,['Estimated time remaining: ' num2str(round(remainingTime)) ' seconds']);
            tic;
            drawnow;
        end;
    end;
end;