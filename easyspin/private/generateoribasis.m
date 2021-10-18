% Generate array of orientational basis functions indices
% 
% Basis = generateoribasis(Basis,basistype)
%
% Input:
%   Basis      structure with basis set information
%    .evenLmax maximum even L
%    .oddLmax  maximum odd L
%    .Kmax     maximum K
%    .Mmax     maximum M
%    .evenK    whether to use only even K (default: false)
%    .jKmin    minimum jK (-1 or 1) (optional; default -1)
%   basistype  string specifying type of basis set
%              'LjKKM','LMK','LMKjK'
%
% Output:
%   Basis.L
%   Basis.K
%   Basis.M
%   Basis.jK
%   
function BasisNew = generateoribasis(Basis,basistype)

if nargin < 2
  basistype = 'LjKKM';
end

Leven = 0:2:Basis.evenLmax;
Lodd  = 1:2:Basis.oddLmax;
Basis.Llist = sort([Leven Lodd]);

if isfield(Basis,'jKmin')
  if Basis.jKmin~=1 && Basis.jKmin~=-1
    error('Basis.jKmin must be either -1 or +1.');
  end
else
  Basis.jKmin = -1;
end

if ~isfield(Basis,'evenK')
  Basis.evenK = false;
end

% Make sure Kmax is even
if Basis.evenK
  Basis.Kmax = floor(Basis.Kmax/2)*2;
end

switch basistype
  case 'LjKKM'
    BasisNew = generatebasis_LjKKM(Basis);
  case 'LMKjK'
    BasisNew = generatebasis_LMKjK(Basis);
  case 'LMK'
    BasisNew = generatebasis_LMK(Basis);
end

return

%-------------------------------------------------------------------------------
% LjKKM: K-symmetrized basis in ordering as in Freed programs:
% 1. increasing L, 2. increasing jK, 3. increasing K, 4. increasing M
%-------------------------------------------------------------------------------
function BasisNew = generatebasis_LjKKM(Basis)

Kmax = Basis.Kmax;
Mmax = Basis.Mmax;

jKmin = Basis.jKmin;
deltaK = Basis.evenK+1;
iBasis = 0;
for L = Basis.Llist
  if mod(L,2)==0, Lparity = +1; else, Lparity = -1; end
  for jK = jKmin:2:1
    Kmx = min(L,Kmax);
    for K = 0:deltaK:Kmx
      if K==0 && (Lparity~=jK), continue; end
      Mmx = min(L,Mmax);
      for M = -Mmx:1:Mmx
        iBasis = iBasis + 1;
        basisList(iBasis,:) = [L jK K M];
      end % M
    end % K
  end % jK
end % L

BasisNew = Basis;
BasisNew.L = basisList(:,1);
BasisNew.jK = basisList(:,2);
BasisNew.K = basisList(:,3);
BasisNew.M = basisList(:,4);

%-------------------------------------------------------------------------------
% LMKjK: K-symmetrized basis in ordering close to LMK
% 1. increasing L, 2. increasing M, 3. increasing K, 4. increasing jK
%-------------------------------------------------------------------------------
function BasisNew = generatebasis_LMKjK(Basis)

Kmax = Basis.Kmax;
Mmax = Basis.Mmax;

jKmin = Basis.jKmin;
deltaK = Basis.evenK+1;
iBasis = 0;
for L = Basis.Llist
  Mmx = min(L,Mmax);
  for M = -Mmx:1:Mmx
    Kmx = min(L,Kmax);
    for K = -Kmx:deltaK:Kmx
      if (K==0), jK = (-1)^L; else, jK= sign(K); end
      if jK<jKmin, continue; end
      iBasis = iBasis + 1;
      basisList(iBasis,:) = [L M abs(K) jK];
    end % K
  end % M
end % L

BasisNew = Basis;
BasisNew.L = basisList(:,1);
BasisNew.M = basisList(:,2);
BasisNew.K = basisList(:,3);
BasisNew.jK = basisList(:,4);

return

%-------------------------------------------------------------------------------
% LMK: LMK basis (not K-symmetrized) in the following ordering: 
% 1. increasing L, 2. increasing M, 3. increasing K
%-------------------------------------------------------------------------------
function BasisNew = generatebasis_LMK(Basis)

Kmax = Basis.Kmax;
Mmax = Basis.Mmax;

deltaK = Basis.evenK+1;
iBasis = 0;
for L = Basis.Llist
  Mmx = min(L,Mmax);
  for M = -Mmx:1:Mmx
    Kmx = min(L,Kmax);
    for K = -Kmx:deltaK:Kmx
      iBasis = iBasis + 1;
      basisList(iBasis,:) = [L M K];
    end % K
  end % M
end % L

BasisNew = Basis;
BasisNew.L = basisList(:,1);
BasisNew.M = basisList(:,2);
BasisNew.K = basisList(:,3);
