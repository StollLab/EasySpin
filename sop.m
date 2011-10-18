% sop  Spin operator matrices
%
%   SpinOp = sop(SpinSystem, Comps)
%   [SpinOp1,SpinOp2,...] = sop(SpinSystem,Comps1,Comps2,...)
%   ... = sop(...,'sparse')
%
%   Spin operator matrix of the spin system
%   SpinSystem in the standard |mS,mI,...>
%   representation.
%
%   If more than one component is given, a matrix
%   is computed for each component.
%
%   Input:
%   - SpinSystem: vector of spin quantum numbers
%     or a spin system specification structure
%   - Comps: string containing 'e','x','y','z','+','-'
%     for each spin, indicating E,Sx,Sy,Sz,S+,S-
%
%   Output:
%   - SpinOp: operator matrix as requested
%
%   Examples:
%     sop([1/2 1],'xy')    % returns SxIy for a S=1/2, I=1 system.
%
%     sop([1/2 1/2],'e+')  % returns SeI+ for a S=I=1/2 system.
%
%     [Sx,Sy,Sz] = sop(1/2,'x','y','z')  % computes three matrices in one go.

function varargout = sop(SpinSystem,varargin)

% UNDOCUMENTED OLD SYNTAX
%  OLD: sop([1/2 1],2,1) equivalent to
%  NEW: sop([1/2 1],'xe')
%
%   OLD SYNTAX  sop(sys,spins,coords)
%     spins: list of indices into sys
%     coords: 1=x, 2=y, 3=z, 4=+, 5=-, all others=e
%
%   DO NOT REMOVE THE HANDLING OF THE OLD SYNTAX!!!!
%   MANY FUNCTIONS RELY ON IT, e.g. zeeman, internal, stev, hfine, resfields

% UNDOCUMENTED THIRD ARGUMENT:
%   - 'sparse' (optional): if specified, SpinOp
%     is returned in sparse form.

% For the spin system J1,J2,J3... the
% order of the spin states of the basis is
%
%  |J1   J2   ...>
%  |J1  J2-1  ...>
%      ...
%  |J1   -J2  ...>
%  |J1-1 J2   ...>
%      ...
%      ...
%  |-J1  -J2  ...>

if (nargin==0), help(mfilename); return; end

if isstruct(SpinSystem)
  SpinVec = spinvec(SpinSystem);
else
  SpinVec = SpinSystem;
end

SparseOut = strcmpi(varargin{end},'sparse');

if SparseOut
  Coords = varargin(1:end-1);
else
  Coords = varargin(1:end);
end

if numel(Coords)==0, error('Not enough input arguments!'); end

OldSyntax = ~ischar(varargin{1});

if OldSyntax,
  %warning(sprintf(['You are using the old syntax of sop (version 1.1 and earlier)!\n'...
  %'It is not guaranteed to work any more! Change to the new syntax!!']));
  
  Spins = varargin{1};
  Coords = varargin{2};
else
  if numel(Coords)>1
    if nargout~=numel(Coords)
      error('Number of input arguments after the first and total number of output arguments do not match.');
    end
    for k=1:numel(Coords)
      if SparseOut
        varargout{k} = sop(SpinVec,Coords{k},'sparse');
      else
        varargout{k} = sop(SpinVec,Coords{k});
      end
    end
    return
  else
    Coords = double(lower(varargin{1}));
    Coords(Coords=='e')=0;
    Coords(Coords=='x')=1;
    Coords(Coords=='y')=2;
    Coords(Coords=='z')=3;
    Coords(Coords=='+')=4;
    Coords(Coords=='-')=5;
    Coords(Coords=='p')=4;
    Coords(Coords=='m')=5;
    if any(Coords>5)
      error('Unrecognized component specification! Must be e, x, y, z, +, or -.');
    end
    Spins = 1:length(Coords);
    if length(Coords)~=length(SpinVec),
      error('Number of spins and number of components do not match!');
    end
  end
end

% identity operator for each spin
Comps = zeros(size(SpinVec));
% and specified components for specified spins
Comps(Spins) = Coords;

% the starting null-spin space is one-dimensional
ia = 1; % 1st (column) index
ja = 1; % 2nd (row) index
sa = 1; % value
na = 1; % dimension

% Run over all spins
for iSpin = 1:length(SpinVec)
  I = SpinVec(iSpin);
  if (I<=0), continue; end

  n = 2*I+1;
  
  % Component switchyard
  %----------------------------------------
  switch Comps(iSpin)
   case 1 % x component
    m = (1:n-1).';
    Dia = 1/2*sqrt(m.*m(end:-1:1));
    ib = [m; m+1];
    jb = [m+1; m];
    sb = [Dia; Dia];
   case 2 % y component
    m = (1:n-1).';
    Dia = -0.5i*sqrt(m.*m(end:-1:1));
    ib = [m; m+1];
    jb = [m+1; m];
    sb = [Dia; -Dia];
   case 3 % z component
    m = (1:n).';
    ib = m;
    jb = m;
    sb = I+1-m;
   case 4 % up shift
    m = (1:n-1).';
    ib = m;
    jb = m+1;
    sb = sqrt(m.*m(end:-1:1));
   case 5 % down shift
    m = (1:n-1).';
    ib = m+1;
    jb = m;
    sb = sqrt(m.*m(end:-1:1));
   otherwise % officially 0, identity
    m = (1:n).';
    ib = m;
    jb = m;
    sb = ones(n,1);
  end
  
  % Kronecker product in sparse form
  % operates only on values and indices.
  %----------------------------------------
  % expansion vectors for new indices and values
  ka = ones(size(sa));
  kb = ones(size(sb));
  
  % new column indices
  t = n*(ia-1).';
  ia = t(kb,:) + ib(:,ka);
  ia = ia(:);
  
  % new row indices
  t = n*(ja-1).';
  ja = t(kb,:) + jb(:,ka);
  ja = ja(:);
  
  % new values
  sa = sb*sa.';
  sa = sa(:);
  
  % new dimension
  na = na*n;
end

% construct sparse matrix
SpinOp = sparse(ia,ja,sa,na,na);

% possibly convert sparse to full matrix, output
if ~SparseOut, SpinOp = full(SpinOp); end

varargout = {SpinOp};

return
