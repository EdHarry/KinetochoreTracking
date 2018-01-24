function success = closeImaris( this )
% IceImarisConnector:  closeImaris (public method)
%
% DESCRIPTION
%
%   This method close the Imaris instance associated to the 
%   IceImarisConnector object (conn) and resets the 
%   mImarisApplication property.
%
% SYNOPSIS
%
%   success = conn.closeImaris( )
%
% INPUT
%
%   None
%
% OUTPUT
% 
%   success : 1 if closing Imaris was successful, 0 otherwise

% LICENSE
%
% ImarisConnector is a simple commodity class that eases communication between 
% Imaris and MATLAB using the Imaris XT interface.
% Copyright (C) 2011  Aaron Ponti
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


% If Imaris is not running, we return here
% We return to signal that Imaris is close, even though we haven't
% closed it
if this.isAlive() == false
    success = 1;
    return
end

try
    this.mImarisApplication.Quit( );
    this.mImarisApplication = [ ];
    success = 1;
catch ex
    disp( ex.message );
    success = 0;
end

end
