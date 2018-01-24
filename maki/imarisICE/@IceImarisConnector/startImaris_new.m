function success = startImaris_new( this, userControl )

% EHarry Feb 2012

%% ORIGINAL HEADER
% % % IceImarisConnector:  startImaris (public method)
% % %
% % % DESCRIPTION
% % %
% % %   This method starts an Imaris instance and stores the ImarisApplication
% % %   ICE object.
% % %
% % %
% % % SYNOPSIS
% % %
% % %   success = conn.startImaris( userControl )
% % %
% % % INPUT
% % %
% % %   userControl : (optional, default = 0 )
% % %                 The optional parameter userControl sets the fate of
% % %                 Imaris when the client is closed: if userControl is
% % %                 true (1), Imaris terminates when the IceImarisConnector
% % %                 object (conn) is deleted. If is it set to false (0),
% % %                 Imaris stays open after the client is closed.
% % %
% % % OUTPUT
% % %
% % %   success : 1 if starting Imaris was successful, 0 otherwise
% %
% % % LICENSE
% % %
% % % ImarisConnector is a simple commodity class that eases communication between
% % % Imaris and MATLAB using the Imaris XT interface.
% % % Copyright (C) 2011  Aaron Ponti
% % %
% % % This program is free software; you can redistribute it and/or
% % % modify it under the terms of the GNU General Public License
% % % as published by the Free Software Foundation; either version 2
% % % of the License, or (at your option) any later version.
% % %
% % % This program is distributed in the hope that it will be useful,
% % % but WITHOUT ANY WARRANTY; without even the implied warranty of
% % % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% % % GNU General Public License for more details.
% % %
% % % You should have received a copy of the GNU General Public License
% % % along with this program; if not, write to the Free Software
% % % Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
% %
% % % Imaris only runs on Windows and Mac OS X
success = 0;
if isempty( strfind( computer, 'PCWIN' ) ) && ...
        isempty( strfind( computer, 'MAC' )  )
    disp( 'IceImarisConnector can only work on Windows and Mac OS X' );
    return
end

% Check and if needed set the optional parameter userControl
if nargin == 1
    userControl = false;
end

if ~ismember( userControl, [ 0 1 ] )
    error( 'userControl must be either 0 or 1' );
end

% Store the userControl
this.mUserControl = userControl;

% If an Imaris instance is open, we close it -- no questions asked
if this.isAlive( ) == 1
    this.closeImaris( );
end

% Now we open a new one
try
%     Launch Imaris
%         [ status, result ] = system( ...
%             [ '"', this.mImarisExePath, '" id', ...
%             num2str( this.mImarisObjectID ), ' &' ] );
%         if status == 1
%             disp( result );
%             success = 0;
%             return
%         end

% Launch Imaris
        [ status, result ] = system(['osascript ~/Documents/MATLAB/maki/imarisICE/startImarisAppleScript.scpt ' num2str( this.mImarisObjectID ) ' &']);
        if status == 1
            disp( result );
            success = 0;
            return
        end


%     imPath = regexprep(this.mImarisExePath,' ','~w');
%     
%     s = ['!/usr/bin/osascript /Users/Ed/Documents/MATLAB/maki/imarisICE/runPythonScript.scpt ',imPath,' ',num2str( this.mImarisObjectID ),' &'];
%     
%     eval(s);
    
    % Add a short pause to let the loading complete
    pause( 2 );
    
    % Try getting the application over a certain time period in case it
    % takes to long for Imaris to be regeistered
    for trial = 1 : 200
        try
            vImaris = this.mImarisLib.GetApplication( this.mImarisObjectID );
        catch ex
            % Silent exception
        end
        if ~isempty( vImaris )
            break
        end
        pause(0.1)
    end
    
    % At this point we should have the application
    if isempty( vImaris )
        disp( 'Could not link to the Imaris application.' );
        success = 0;
        return
    end
    
    % We can store the application
    this.mImarisApplication = vImaris;
    
    % We set success to 1
    success = 1;
    
catch ex
    
    disp( ex.message );
    success = 0;
    
end

end
