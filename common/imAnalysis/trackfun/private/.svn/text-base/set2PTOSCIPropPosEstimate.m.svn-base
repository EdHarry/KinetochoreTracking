function set2PTOSCIPropPosEstimate(pos)
% updates the current position and parameter vector of the 2PTOSCI 
% coordinate propagator
global prop__;

% store previous position
if(~isempty(prop__.currentPos))
   prevPos = prop__.currentPos;
else
   prevPos = [];
end;

% update position
prop__.currentPos = pos;

% compute new azimuth and excursion
if(prevPos)
   crds = [prevPos;prop__.currentPos];
   dd = diff(crds);
   newExc = norm(dd);
   newAzi = atan2(dd(2),dd(1));
   if(newAzi < 0)
      newAzi = 2*pi + newAzi;
   end;
   if(isempty(prop__.params))
      prop__.params(1) = newAzi;
      prop__.params(2) = newExc;
      prop__.params(3) = 1;
   else
      % update the parameters as the mean of old and new data
      prop__.params(1) = (prop__.params(3)*prop__.params(1) + newAzi)/...
         (prop__.params(3) + 1);
      prop__.params(2) = (prop__.params(3)*prop__.params(2) + newExc)/...
         (prop__.params(3) + 1);
      prop__.params(3) = prop__.params(3) + 1;
   end;
   
   % flip the azimuth because of the oscillation
   prop__.params(1) = mod(prop__.params(1)+pi,2*pi);
   prop__.params;
end;

   
