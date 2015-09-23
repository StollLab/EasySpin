% sxliou takes the Hilbert space spin operator Sx as input and returns as
% output the vector representation of the Sx operator in Liouville space

function SxL = startvec(basisList,SxH)


Lmax = max(basisList(:,1));
SpaceSize = sum((2*(0:Lmax)+1).^2);

P_Spatial = zeros(size(basisList,1),1); % length takes longest dim!!!!
%P_Spatial(1) = sqrt(1/(8*pi^2))
P_Spatial(1) = 1; % assumes L=M=K=0 is the first basis function!

% I think this is incorrect for more than one electron spin
%{
nSpins = numel(SxH);
for iS = 1:nSpins
  
  % normalize the starting vector
  
  %SxVector = SxH{iS}(:)/norm(SxH{iS}(:));
  SxVector = SxH{iS}(:);
  SxL = kron(P_Spatial,SxVector);
  %SxL{iS} = kron(SxH{iS}(:),P_Spatial);
  
end

for iS = 1:nSpins
  
  % normalize the starting vector
  SxVector = SxH{iS}(:)/norm(SxH{iS}(:));
  SxL = kron(P_Spatial,SxVector);
  %SxL{iS} = kron(SxH{iS}(:),P_Spatial);
  
end

%}

SxVector = SxH(:);
SxVector = SxVector/norm(SxVector);
SxL = kron(P_Spatial,SxVector);
SxL = sparse(SxL);

return