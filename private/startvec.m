% startvec takes the Hilbert space spin operator Sx as input and returns as
% output the vector representation of the Sx operator in Liouville space

function SxL = startvec(basisList,SxH)

% assumes L=M=K=0 is the first spatial basis function!
SxVector = SxH(:)/norm(SxH(:));
nSpin = numel(SxVector);
nBasis = nSpin*size(basisList,1);
SxL = sparse(1:nSpin,1,SxVector,nBasis,nBasis);

return
