function cMap = getCustomCMap(id)
%getCustomCMap: Get customized color map.
%
% SYNOPSIS: cMap = getCustomCMap(id)
%
% INPUT:
%    id: Color map id:
%        1: hot color: Green to yellow to red to dark red;
%        2: cold color: Dark blue to blue to cyan;
%        3: mixture of hot and cold (1 & 2) to distinguish between two score range.
%        4: Blue to green to yellow to red.
%        5: 'jet' map but without dark red.
%        6: 'jet' map but with shortened cold color range to increase the
%           contrast from cold to hot.
%        7: First half: blue to green; Second half: Yellow to red.

%Hot or positive color.
%Green to yellow.
hotMap = [linspace(0.5,1,32).' ones(32,1) zeros(32,1)];
%Yellow to red.
hotMap = [hotMap; [ones(32,1) linspace(1,0,32).' zeros(32,1)]];
%Red to dark red.
hotMap = [hotMap; [linspace(1,0.5,16).' zeros(16,1) zeros(16,1)]];

%Cold or negative color.
%Dark blue to blue.
coldMap = [zeros(16,1) zeros(16,1) linspace(0.5,1,16).'];
%Blue to cyan.
coldMap = [coldMap; [zeros(64,1) linspace(0,1,64).' ones(64,1)]];

%Jet map but without dark red.
%blue to green.
bgyrMap = [zeros(32,1) linspace(0,1,32).' linspace(1,0,32).'];
%Green to yellow.
bgyrMap = [bgyrMap; linspace(0.5,1,32).' ones(32,1) zeros(32,1)];
%Yellow to red.
bgyrMap = [bgyrMap; [ones(32,1) linspace(1,0,32).' zeros(32,1)]];

%Jet map.
jetMap = colormap(jet(148));

switch id
   case 1
      cMap = hotMap;
   case 2
      cMap = coldMap;
   case 3 
      %Mixture of case 1 and 2. Ideal for distinguishing two score ranges.
      %For example, poly/depoly or protrusion/retraction.
      cMap = [coldMap; hotMap];
   case 4
      cMap = bgyrMap;
   case 5
      cMap = jetMap(1:128,:);
   case 6
      cMap = zeros(118,3);
      cMap(27:118,:) = jetMap(57:148,:);
      cMap(1:5,:) = [zeros(5,1) zeros(5,1) linspace(0.5135,1,5).'];
      cMap(6:26,:) = [zeros(21,1) linspace(0.0270,1,21).' ones(21,1)];
   case 7
      cMap = [zeros(64,1) linspace(0,1,64).' linspace(1,0,64).'];
      cMap = [cMap; [ones(64,1) linspace(1,0,64).' zeros(64,1)]];
end
