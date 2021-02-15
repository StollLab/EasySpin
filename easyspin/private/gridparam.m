% gridparam  Get grid parameter for a given point group symmetry
%
%   maxPhi = gridparam(PointGroup)
%   [maxPhi,closedPhi] = gridparam(PointGroup)
%   [maxPhi,closedPhi,nOctants] = gridparam(PointGroup)
%
% Inputs:
%   PointGroup: the Schoenfliess name of a centrosymmetric group or C1
%
% Outputs:
%   maxPhi:      the maximum phi of the unique spherical surface
%   closedPhi:   open or closed phi interval
%   nOctants:    how many octants sphgrid uses for this group

% Group  maxPhi   nOctants   closedPhi
% C1     2*pi     8          false
% Ci     2*pi     4          false
% C2h    pi       2          false
% S6     2*pi/3   2          false   (could also work with 1 octant)
% C4h    pi/2     1          false
% C6h    pi/3     1          false
% D2h    pi/2     1          true
% Th     pi/2     1          true
% D3d    pi/3     1          true
% D4h    pi/4     1          true
% Oh     pi/4     1          true
% D6h    pi/6     1          true
% Dinfh  0        0          true
% O3     0        -1         true

function varargout = gridparam(PointGroup)

maxPhi = [2 2 1 2/3 1/2 1/3 1/2 1/2 1/3 1/4 1/4 1/6 0 0]*pi;
nOctants = [8 4 2 2 1 1 1 1 1 1 1 1 0 -1];
closedPhi = [false false false false false false true true true true ...
           true true true true];
Groups = {'C1','Ci','C2h','S6','C4h','C6h','D2h','Th',...
          'D3d','D4h','Oh','D6h','Dinfh','O3'};        
        
idx = find(strcmp(PointGroup,Groups));
if isempty(idx)
  error('Unsupported symmetry group %s!',PointGroup);
else
  varargout = {maxPhi(idx),closedPhi(idx),nOctants(idx)};
  varargout = varargout(1:nargout);
end
