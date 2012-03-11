% symparam  Get symmetry group parameters [EasySpin/private]
%
%   maxPhi = symparam(GroupName)
%   [maxPhi,openPhi] = symparam(GroupName)
%   [maxPhi,openPhi,nOctants] = symparam(GroupName)
%
%   GroupName: The Schoenfliess name of a centrosymmetric group
%   maxPhi: the maximum phi of the unique spherical surface
%   openPhi: open or closed phi interval?
%   nOctants: how many octants sphgrid uses for this group

function varargout = symparam(SymGroup)

maxPhi = [2 2 1 2/3 1/2 1/3 1/2 1/2 1/3 1/4 1/4 1/6 0 0]*pi;
nOctants = [8 4 2 2 1 1 1 1 1 1 1 1 0 -1];
openPhi = [1 1 1 1 1 1 0 0 0 0 0 0 0 0];
Groups = {'C1','Ci','C2h','S6','C4h','C6h','D2h','Th',...
          'D3d','D4h','Oh','D6h','Dinfh','O3'};

idx = strcmp(SymGroup,Groups);
if isempty(idx)
  error('Unsupported symmetry group %s!',SymGroup);
else
  varargout = {maxPhi(idx),openPhi(idx),nOctants(idx)};
  varargout = varargout(1:nargout);
end
