function createBlankGui(filename)
% createBlankGui is a work-around for MATLAB's bug 33257 (MATLAB v6.5)
%
% GUIDE fails to open one of its UI templates and crashes in one of several ways (see
%    http://www.mathworks.com/support/solutions/data/33257.shtml for more details).
%
% SYNOPSIS      createBlankGui(filename)
%
% INPUT         filename    : (optional) File name (with full path) of the GUI to be created. 
%                             The extension .fig is optional (if missing is optionally added). 
%
%                             If createBlankGUI is called with no arguments, the user will be
%                             asked to specify the file name through a save dialog.
%
%                             GUIDE will be automatically launched to further edit the new GUI.
%
% OUTPUT        none
%
% DEPENDENCES   createBlankGui uses { }
%               createBlankGui is used by { }
%
% Aaron Ponti, June 7th, 2004

% Check input
if nargin==0
    
    % Ask the user to specify a file name for the GUI
    [filename,path]=uiputfile('','Specify path and name of the GUI you want to create');
    
    if isequal(path,0) | isequal(filename,0)
        
        % The user pressed 'Cancel'
        msgbox('Aborted by the user.','Info','warn');
        return
        
    end
    
end

% Create a figure and save it 
h=figure;
try
    
    % Try and save the figure
    hgsave(h,[path,filesep,filename]);
    
    % Close the figure
    close(h);
    
catch

    % Close the figure
    close(h);
    
    % Return an error and leave
    errordlg('The figure could not be saved.','Error');
    return
    
end
 
% Open it with GUIDE
try
    
    % Try and open the figure
    guide([path,filesep,filename]);
    
catch
    
    % Return an error and leave
    errordlg('The created figure could not be loaded by GUIDE.','Error');
    return

end
