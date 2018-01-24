% IceImarisConnector constructor
%
% DESCRIPTION
%
%   IceImarisConnector is a simple commodity class that eases communication
%   between Imaris and MATLAB using the Imaris XT interface.
%
% SYNOPSIS
%
%   conn = IceImarisConnector( vImarisApplication )
%
% INPUT
%
%   vImarisApplication : (optional) if omitted, an IceImarisConnector
%                        object is created that is not connected to any 
%                        Imaris instance.
%                        Imaris can then be started (and connected) using
%                        the startImaris( ) method, i.e.
%
%                            conn.startImaris( )
%
%                        Alternatively, vImarisApplication can be an
%                        Imaris Application ID (as provided by Imaris)
%                        or directly an Imaris Application ICE object.
%
% REMARK
%
%   The Imaris Application ICE object is stored in the read-only property
%   mImarisApplication. The mImarisApplication property gives access to
%   the entire Imaris ICE methods. Example:
%
%   conn.mImarisApplication.GetSurpassSelection( )
%
%   returns the currently selected object in the Imaris surpass scene.
%
% OUTPUT
%
%   conn    :  an object of class IceImarisConnector

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

classdef IceImarisConnector < handle
        
    % Properties
    properties ( SetAccess = private, GetAccess = public )
        mImarisApplication = [];
    end
    
    properties ( Access = private )
        mImarisLib;
        mImarisExePath;
        mImarisLibPath;
        mImarisObjectID;
        mUserControl;
    end
    
    % Methods
    methods
        
        % Constructor
        function this = IceImarisConnector( vImarisApplication )

            % First, we prepare everything we need
            % Store the Imaris and ImarisLib path
            [ success, errorMessage ] = findImaris( this );
            if ~success
                error( errorMessage );
            end
            
            % Add the ImarisLib.jar package to the java class path
            % (if not there yet)
            if all( cellfun( ...
                    @isempty, strfind( javaclasspath, 'ImarisLib.jar' ) ) )
                javaaddpath( this.mImarisLibPath );
            end
            
            % Create and store an ImarisLib instance
            this.mImarisLib = ImarisLib( );

            % Assign a random id
            this.mImarisObjectID = randi( 100000 );

            % Now we check the (optional) input parameter.
            % We have three cases. If called without parameters, we
            % just create an IceImarisConnector object that does nothing.
            % If we get one input parameter, we have to distinguish between
            % two cases: We either get an ID (as is the case when the 
            % function is launched from Imaris) and thus we query for 
            % the application, or directly the object: in this case we 
            % just assign it to the mImarisApplication property
            
            if nargin == 0
                
                % We already did everything
                return
                
            elseif nargin == 1
                
                if isa( vImarisApplication, 'IceImarisConnector' )
                    
                    % If the input parameter is an IceImarisConnector
                    % object we return the reference. This way, an
                    % XTension class can take a reference to an 
                    % IceImarisConnector object as input parameter
                    this = vImarisApplication;
                    
                elseif isa( vImarisApplication, 'Imaris.IApplicationPrxHelper')
                    
                    % This is an Imaris application object - se store it
                    this.mImarisApplication = vImarisApplication;
                    
                elseif isscalar( vImarisApplication )

                    % Check if the application is registered
                    server = this.mImarisLib.GetServer;
                    nApps = server.GetNumberOfObjects( );
                    if nApps == 0
                        error( 'There are no registered Imaris applications.' );
                    end
                    
                    if server.GetObjectID( vImarisApplication ) == -1
                        error( 'Invalid Imaris application ID.' );
                    end

                    this.mImarisApplication = ...
                        this.mImarisLib.GetApplication( vImarisApplication );
                
                else
                
                    error( [ 'The passed object is not a Imaris ', ...
                        'Application ID.' ] );
   
                end
            
            else
                
                error( 'Wrong number of input arguments!' );
                
            end
            
        end
       
        % Desctructor
        function delete( this )
        
            if this.mUserControl == 1
                if ~isempty( this.mImarisApplication )
                    this.closeImaris( );
                end
            end

        end
        
    end
    
    % External public and static methods
    % =====================================================================
    
    methods ( Access = public )
        
        % display
        display( this )        

        % getAllSurpassChildren
        children = getAllSurpassChildren( this, recursive )

        % getDataVolume
        stack = getDataVolume( this, channel, timepoint, iDataset )
        
        % getExtends
        extends = getExtends( this )

        % getImarisVersionAsInteger
        version = getImarisVersionAsInteger( this )
        
        % getMatlabDatatype
        type = getMatlabDatatype( this )
        
        % getSizes
        sizes = getSizes( this )

        % getVoxelSizes
        voxelSizes = getVoxelSizes( this )
        
        % isAlive
        alive = isAlive( this )
        
        % mapPositionsUnitsToVoxels
        varargout = mapPositionsUnitsToVoxels( this, varargin )
        
        % mapPositionsVoxelsToUnits
        varargout = mapPositionsVoxelsToUnits( this, varargin )
        
        % startImaris
        success = startImaris( this, userControl )
        
        % getSurpassCameraRotationMatrix
        [ R, isI ] = getSurpassCameraRotationMatrix( this )

        % setDataVolume
        setDataVolume( this, stack, channel, timepoint )

    end

    methods ( Access = public, Static = true )
        
        % getVersion
        function version = getVersion( )
            version = '0.0.1';
        end
        
    end
    
    methods ( Access = private )
        
        % findImaris
        [ imarisPath, errorMessage ] = findImaris( this )
            
    end
    
end