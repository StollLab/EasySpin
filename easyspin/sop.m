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
%     SxIy = sop([1/2 1],'xy')    % returns SxIy for a S=1/2, I=1 system.
%
%     SeIp = sop([1/2 1/2],'e+')  % returns SeI+ for a S=I=1/2 system.
%
%     [Sx,Sy,Sz] = sop(1/2,'x','y','z')  % computes three matrices in one go.

function varargout = sop(SpinSystem,varargin)

% UNDOCUMENTED OLD SYNTAX
%  OLD: sop([1/2 1],2,1) equivalent to
%  NEW: sop([1/2 1],'xe')
%
%   OLD SYNTAX  sop(sys,spins,comps)
%     spins: list of indices into sys
%     comps: 1=x, 2=y, 3=z, 4=+, 5=-, all others=e
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
nSpins = numel(SpinVec);

SparseOutput = strcmpi(varargin{end},'sparse');

if SparseOutput
  Coords = varargin(1:end-1);
else
  Coords = varargin(1:end);
end

if numel(Coords)==0
  error('Not enough input arguments!');
end

OldSyntax = ~ischar(varargin{1});

if OldSyntax,
  
  Spins = varargin{1};
  Coords = varargin{2};

else
  
  if nargout~=numel(Coords)
    error('Number of output arguments (%d) and number of requested operator matrices (%d) do not match.',nargout,numel(Coords));
  end
  
  if numel(Coords)>1
    for k=1:numel(Coords)
      if SparseOutput
        varargout{k} = sop(SpinVec,Coords{k},'sparse');
      else
        varargout{k} = sop(SpinVec,Coords{k});
      end
    end
    return
  else
    Coords = lower(varargin{1});
    
    % Check syntax of component string
    syntax1 = sprintf('^[exyz+\\-pmab]{%d}$',nSpins);
    isSyntax1 = ~isempty(regexp(Coords,syntax1,'match'));
    if ~isSyntax1
      syntax2 = sprintf('^([exyz+\\-pmab]\\d+){1,%d}$',nSpins);
      isSyntax2 = ~isempty(regexp(Coords,syntax2,'match'));
      if ~isSyntax2
        error('Could not determine what ''%s'' is for the given spin system.',Coords);
      end
    end
   
    if isSyntax1
      Coords = double(Coords);
      Spins = 1:nSpins;
    elseif isSyntax2
      % Get tokens from pattern
      pattern = '([exyz+\-pmab]\d+)+?';
      tokens = char(regexp(Coords,pattern,'match'));
      % Get list of spin indices
      Spins = str2num(tokens(:,2:end)).'; %#ok<ST2NM>
      if any(Spins>nSpins)
        error('sop: The spin system only contains %d spins, but you requested ''%s''.',nSpins,Coords);
      end
      if numel(unique(Spins))<numel(Spins)
        error('sop: Repeated spin in %s',Coords');
      end
      % Get list of requested components
      Coords = double(tokens(:,1).');
    end

    if (any(SpinVec(Coords=='a')~=1/2)) || ...
       (any(SpinVec(Coords=='b')~=1/2))
      error('''a'' and ''b'' work only for spin-1/2.');
    end
        
  end
end

% identity operator for each spin
Comps = repmat('e',1,nSpins);
% and specified components for specified spins
Comps(Spins) = Coords;

% the starting null-spin space is one-dimensional
ia = 1; % 1st (column) index
ja = 1; % 2nd (row) index
sa = 1; % value
na = 1; % dimension

% Run over all spins
for iSpin = 1:nSpins
  I = SpinVec(iSpin);
  if (I<=0), continue; end

  n = 2*I+1;
  
  % Component switchyard
  %----------------------------------------
  switch Comps(iSpin)
    case {'x',1} % x component
      m = (1:n-1).';
      Dia = 1/2*sqrt(m.*m(end:-1:1));
      ib = [m; m+1];
      jb = [m+1; m];
      sb = [Dia; Dia];
    case {'y',2} % y component
      m = (1:n-1).';
      Dia = -0.5i*sqrt(m.*m(end:-1:1));
      ib = [m; m+1];
      jb = [m+1; m];
      sb = [Dia; -Dia];
    case {'z',3} % z component
      m = (1:n).';
      ib = m;
      jb = m;
      sb = I+1-m;
    case {'+','p',4} % up shift
      m = (1:n-1).';
      ib = m;
      jb = m+1;
      sb = sqrt(m.*m(end:-1:1));
    case {'-','m',5} % down shift
      m = (1:n-1).';
      ib = m+1;
      jb = m;
      sb = sqrt(m.*m(end:-1:1));
    case 'a' % alpha, for spin-1/2 only
      ib = 1;
      jb = 1;
      sb = 1;
    case 'b' % beta, for spin-1/2 only
      ib = 2;
      jb = 2;
      sb = 1;
    case {'e',0} % identity
      m = (1:n).';
      ib = m;
      jb = m;
      sb = ones(n,1);
    otherwise
      error('Unknown operator specification.');
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
OperatorMatrix = sparse(ia,ja,sa,na,na);

% possibly convert sparse to full matrix, output
if ~SparseOutput
  OperatorMatrix = full(OperatorMatrix);
end

varargout = {OperatorMatrix};

return
