function colormap = isomorphicColormap(color, cLength, test)
%ISOMORPHICCOLORMAP generates isomorphic colormaps
%
% SYNOPSIS: colormap = isomorphicColormap(color,cLength, test)
%
% INPUT color (opt):  'red' (or 'r'), 'green' (or 'g'), 'blue' (or 'b')
%                         for colormaps going in luminosity from
%                         black->COLOR->white
%                     'bw','gw','rw' for blue, green, red changing in
%                         luminosity from b/g/r to white
%                     'b/y' for saturation from blue->grey->yellow
%                     'g/r' for saturation from green->grey->red
%                     'b/o' for saturation from blue->grey->orange
%                     'bko' sat/lum blue-black-orange
%                     'bwo' sat/lum blue-white-orange
%                     'kbwok' lum-sat/lum-lum black-blue-white-orange-black
%                     'kbok' as 'kbwok', but with a sat-jump at 0 to
%                           segment negative vs positive
%
%       cLength (opt): Number of entries of the colormap. Default: 64
%       test (opt) : if 1, colormapTest is being run with the colormap.
%                    Default: 0
% OUTPUT colormap: colormap (64x3 array of RGB values)
%
%
% REMARKS see
%           http://www.research.ibm.com/people/l/lloydt/color/color.htm
%         and
%           http://www.research.ibm.com/dx/proceedings/pravda/index.htm
%         for details
%
% created with MATLAB ver.: 7.3.0.267 (R2006b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 02-Apr-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% defaults
def_color = 'red';
def_test = false;
def_cLength = 64;

% test input
if nargin < 1 || isempty(color)
    color = def_color;
else
    % no need to check here - we will do that with the switch command below
end
if nargin < 2 || isempty(cLength)
    cLength = def_cLength;
end
if nargin < 3 || isempty(test)
    test = def_test;
end


% select colormap
switch color
    case {'red','r'}
        % change color, saturation a tiny bit, luminosity lots
        hue = [linspace(0.95,1,cLength/2)';linspace(0,0.05,cLength/2)'];
        sat = linspace(0.7,0.8,cLength)';
        lum = linspace(0.1,1,cLength)';

        cmap = hsl2rgb([hue,sat,lum]);

    case {'rw'}
        % same as 'r', but going from white to color
        hue = [linspace(0.95,1,cLength/2)';linspace(0,0.05,cLength/2)'];
        sat = repmat(1,cLength,1);
        lum = linspace(1,0.5,cLength)';

        cmap = hsl2rgb([hue,sat,lum]);

    case {'green','g'}
        % change color, saturation a tiny bit, luminosity lots
        hue = linspace(0.3,0.4,cLength)';
        sat = linspace(0.7,0.8,cLength)';
        lum = linspace(0.1,1,cLength)';

        cmap = hsl2rgb([hue,sat,lum]);

    case {'gw'}
        % same as 'g', but going from white to color
        hue = linspace(0.3,0.4,cLength)';
        sat = repmat(1,cLength,1);
        lum = linspace(1,0.5,cLength)';

        cmap = hsl2rgb([hue,sat,lum]);

    case {'blue','b'}
        % change color, saturation a tiny bit, luminosity lots
        hue = linspace(0.5,0.6,cLength)';
        sat = linspace(0.7,0.8,cLength)';
        lum = linspace(0.1,1,cLength)';

        cmap = hsl2rgb([hue,sat,lum]);

    case {'bw'}
        % same as 'g', but going from white to color
        hue = linspace(0.5,0.6,cLength)';
        sat = repmat(1,cLength,1);
        lum = linspace(1,0.5,cLength)';

        cmap = hsl2rgb([hue,sat,lum]);

    case 'b/y'
        % color: blue->yellow, saturation 1->0->1, luminosity 0.5 (=max
        % color)
        hue = repeatEntries([0.66;0.16],cLength/2);
        sat = [linspace(1,0,cLength/2)';linspace(0,1,cLength/2)'];
        lum = ones(cLength,1)*0.50;

        cmap = hsl2rgb([hue,sat,lum]);

    case 'g/r'
        % color: green->red, saturation 1->0->1, luminosity 0.5 (=max
        % color)
        hue = repeatEntries([0.33;1],cLength/2);
        sat = [linspace(0.8,0,cLength/2)';linspace(0,0.8,cLength/2)'];
        lum = ones(cLength,1)*0.50;

        cmap = hsl2rgb([hue,sat,lum]);
    case 'b/o'
        hue = repeatEntries([0.58;0.08],cLength/2);
        sat = [linspace(1,0,cLength/2)';linspace(0,1,cLength/2)'];
        lum = ones(cLength,1)*0.50;

        cmap = hsl2rgb([hue,sat,lum]);
        
            case 'bwo'
        hue = repeatEntries([0.58;0.08],cLength/2);
        sat = [linspace(1,0,cLength/2)';linspace(0,1,cLength/2)'];
        lum = [linspace(0.5,0.9,cLength/2)';linspace(0.9,0.5,cLength/2)'];

        cmap = hsl2rgb([hue,sat,lum]);
         case 'bko'
        hue = repeatEntries([0.58;0.08],cLength/2);
        sat = [linspace(1,0,cLength/2)';linspace(0,1,cLength/2)'];
        lum = [linspace(0.5,0.1,cLength/2)';linspace(0.1,0.5,cLength/2)'];

        cmap = hsl2rgb([hue,sat,lum]);
    case 'kbwok'
        hue = repeatEntries([0.58;0.08],cLength/2);
        sat = [ones(cLength/4,1);linspace(1,0,cLength/4)';linspace(0,1,cLength/4)';ones(cLength/4,1)];
        lum = [linspace(0.1,0.5,cLength/4)';linspace(0.5,0.9,cLength/4)';linspace(0.9,0.5,cLength/4)';linspace(0.5,0.1,cLength/4)'];

        cmap = hsl2rgb([hue,sat,lum]);
        case 'kbok'
        hue = repeatEntries([0.58;0.08],cLength/2);
        sat = [ones(cLength/4,1);linspace(1,0.5,cLength/4)';linspace(0.5,1,cLength/4)';ones(cLength/4,1)];
        lum = [linspace(0.1,0.5,cLength/4)';linspace(0.5,0.75,cLength/4)';linspace(0.75,0.5,cLength/4)';linspace(0.5,0.1,cLength/4)'];

        cmap = hsl2rgb([hue,sat,lum]);
%         case 'go'
%         hue = linspace(0.30,0.1,cLength)';
%         sat = [linspace(0.5,1,cLength)'];
%         %lum = ones(cLength,1)*0.50;
%         lum = [linspace(0.4,0.6,cLength)'];
% 
%         cmap = hsl2rgb([hue,sat,lum]);
    otherwise
        error('color %s not implemented yet',color)
end

% check test
if test
    colormapTest(cmap,color);
end

% check for output argument
if nargout > 0
    colormap = cmap;
end