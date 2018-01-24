function map = QueryCMap(imT)
%IMADJUSTTOOL/QUERYCMAP gets colormap
%
% SYNOPSIS QueryCMap(imT)
%
% INPUT  imT: an ImAdjustTool object
% 
% OUTPUT map : colormap that is associated with current gamma setting;
%
% c: 29/8/00	 gD

map = imT.newCMap;