% stev  Extended Stevens spin operator matrices
%
%   Op = stev(S,k,q)
%   Op = stev(S,k,q,iSpin)
%   Op = stev(S,k,q,iSpin,'sparse')
%
%   Constructs extended Stevens operator matrices for
%   0<=k<=2*S and -k<=q<=k for the spin S.
%
%   If S is a vector representing the spins of a
%   spin system, Op is computed for the spin number
%   iSpin (e.g. the second if iSpin==2) in the state
%   space of the full spin system. If iSpin is omitted,
%   iSpin is set to 1.
%
%   All k values from 0 to 12 are supported. The
%   most common ones are 2, 4 and 6.
%
%   The extended Stevens operators are tesseral (as
%   opposed to spherical) tensor operators and are
%   therefore all Hermitian.
%
%   Input:
%   - S: spin quantum number, or vector thereof
%   - k,q: indices specifying O_k^q
%   - iSpin: index of the spin in the spin vector for
%       which the operator matrix should be computed
%
%   Output:
%   - Op: extended Stevens operator matrix
%
%   Examples:
%    To obtain O_4^2 for a spin 5/2, type
%       stev(5/2,4,2)
%    To obtain O_6^5 for the second spin in a spin
%    system with two spins-3/2, type
%       stev([3/2 3/2],6,5)

% Abbreviations:
%  ESO   extended Stevens operator
%  STO   spherical tensor operator
% Computation of matrix elements based on
%   basic equations for spherical tensor components.
%   See Ryabov, J.Magn.Reson. 140, 141-145 (1999),
%   eqns. 1, 2, 21, 22
% References
%  I.D.Ryabov, J.Magn.Reson. 140, 141-145 (1999)
%  C. Rudowicz, C.Y.Chung, J.Phys.:Condens.Matter 16, 1-23 (2004)
%
% Consistent with
%   Altshuler/Kozyrev, Electron Paramagnetic Resonance
%   in Compounds of Transition Elements, Wiley, 2nd edn., (1974)
%   Appendix V. p. 512
%   Table 16, p. 863 in Abragam/Bleaney
%   Electron Paramagnetic Resonance of Transition Ions, Dover (1986)
% See also
%   St. Stoll, PhD thesis, ETH Zurich, 2003
%   C.Rudowicz, Magn.Reson.Rev. 13, 1-87 (1987)

function Op = stev(Spins,k,q,iSpin,Sparse)

if (nargin==0), help(mfilename); return; end

if (nargin<3) || (nargin>5), error('Wrong number of input arguments!'); end

if (nargin<4)
  iSpin = 1;
end
if (nargin<5)
  Sparse = '';
end

if isstruct(Spins)
  Spins = spinvec(Spins);
end

if ischar(Sparse)
  useSparseMatrices = strcmp(Sparse,'sparse');
else
  if isempty(Sparse)
    useSparseMatrices = false;
  elseif ~islogical(Sparse)
    error('Input argument Sparse must be either true/false or ''sparse''!');
  else
    useSparseMatrices = Sparse;
  end
end


% Initialization of prefactors
%-------------------------------------------------
% Computed in Mathematica using expression from
% Ryabov, J.Magn.Reson. 140, 141-145 (1999)
persistent F;
if isempty(F)
  F(13,1:13) = [1916006400 958003200 958003200 31933440 3991680 1995840 31680 15840 1584 264 24 12 1];
  F(12,1:12) = [319334400 79833600 79833600 13305600 2661120 23760 7920 1320 1320 22 22 1];
  F(11,1:11) = [14515200 7257600 1209600 604800 86400 2880 360 180 20 10 1];
  F(10,1:10) = [1451520 725700 725700 60480 60480 864 288 18 18 1];
  F( 9,1:9)  = [80640 40320 40320 6720 672 336 16 8 1];
  F( 8,1:8)  = [40320 5040 1680 168 168 14 14 1];
  F( 7,1:7)  = [2880 1440 360 60 12 6 1];
  F( 6,1:6)  = [480 240 240 10 10 1];
  F( 5,1:5)  = [48 24 8 4 1];
  F( 4,1:4)  = [24 6 6 1];
  F( 3,1:3)  = [4 2 1];
  F( 2,1:2)  = [2 1];
  F( 1,1)    = 1;
end

% Checks on input parameters S, k and q
%-------------------------------------------------
if (iSpin<0) || (iSpin>numel(Spins))
  error('iSpin = %d is out of range. It should be between 1 and %d',iSpin,numel(Spins));
end

S = Spins(iSpin);

if any(mod(2*Spins,1)) || any(~isreal(Spins))
  error('S must contain positive multiples of 1/2.');
end

kmax = size(F,1)-1;
if numel(k)~=1 || numel(q)~=1
  error('k and q must be single numbers.');
end
if mod(k,1) || (k<0) || (k>kmax) || ~isreal(k)
  error('k too large. Maximum supported k is %d.',kmax);
end
if (k>2*S)
  error('k must not be larger than 2*S (%d with S=%g).',2*S,S);
end
if mod(q,1) || abs(q)>k || ~isreal(q)
  error('q must be an integer between -k and k (%d and %d).',-k,k);
end

% Computation of operator matrix
%-------------------------------------------------
% Compute component of STO using Racah's commutation
% rule, but leaving out scaling and normalization
% constants. This is possible since they are divided
% out again by the Stevens prefactors c (see below).
% (Ryabov, Eq.[1])
Jp = sop(Spins,iSpin,4,'sparse');
Jm = sop(Spins,iSpin,5,'sparse');
T = Jp^k;  % T(k,k)
for qq = k-1:-1:abs(q)
  T = Jm*T - T*Jm;
end

% Linear combination coefficient for the ESO.
% alpha as defined by Ryabov (1999).
alpha = 1;
if mod(q,2) && ~mod(k,2)  % for even k and odd q
  alpha = 1/2;
end
c = alpha/F(k+1,abs(q)+1); % already without N_kq

% The sign of the original normalization constant has
% to be retained.
c = (-1)^(k-q)*c;

% Simple construction formulae for ESOs: construction
% of cosine and sine tesseral operators from STOs.
% (Ryabov, eq.21)
if (q>=0)
  Op = c/2 *(T + T'); % cosine tesseral operator
else
  Op = c/2i*(T - T'); % sine tesseral operator
end

if ~useSparseMatrices
  Op = full(Op);
end

return
