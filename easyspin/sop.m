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
  OperatorSpec = varargin(1:end-1);
else
  OperatorSpec = varargin(1:end);
end

if numel(OperatorSpec)==0
  error('Not enough input arguments!');
end

OldSyntax = ~ischar(varargin{1});

if OldSyntax
  
  Spins = varargin{1};
  Coords = varargin{2};
  
else
  
  if nargout~=numel(OperatorSpec) && nargout>0
    error('Number of output arguments (%d) and number of requested operator matrices (%d) do not match.',nargout,numel(OperatorSpec));
  end
  
  if numel(OperatorSpec)>1
    
    for k = 1:numel(OperatorSpec)
      if SparseOutput
        varargout{k} = sop(SpinVec,OperatorSpec{k},'sparse');
      else
        varargout{k} = sop(SpinVec,OperatorSpec{k});
      end
    end
    return
    
  else
    
    OperatorSpec = lower(varargin{1});
    
    syntax1 = sprintf('^[exyzab+-]{%d}$',nSpins);
    isSyntax1 = ~isempty(regexp(OperatorSpec,syntax1,'match'));
    if ~isSyntax1
      syntax2 = sprintf('^([exyzab+-]\\d+){1,%d}$',nSpins);
      isSyntax2 = ~isempty(regexp(OperatorSpec,syntax2,'match'));
      if ~isSyntax2
        OperatorSpec = regexp(OperatorSpec,',','split');
        syntax3 = '^[exyzab+-](\(\d+(\|\d+)?\))?(\d+)?$';
        matched = regexp(OperatorSpec,syntax3,'match');
        isSyntax3 = ~any(cellfun(@isempty, matched));
        if ~isSyntax3
          error('Could not determine what ''%s'' is for the given spin system.',lower(varargin{1}));
        end
      end
    end
          
    if isSyntax1
      
      Coords = double(OperatorSpec);
      Spins = 1:nSpins;
      
    elseif isSyntax2
      
      % Get tokens from pattern
      pattern = '([exyzab+-]\d+)+?';
      tokens = char(regexp(OperatorSpec,pattern,'match'));
      % Get list of spin indices
      Spins = str2num(tokens(:,2:end)).'; %#ok<ST2NM>
      if any(Spins>nSpins)
        error('sop: The spin system only contains %d spins, but you requested ''%s''.',nSpins,OperatorSpec);
      end
      if numel(unique(Spins))<numel(Spins)
        error('sop: Repeated spin in %s',OperatorSpec);
      end
      % Get list of requested components
      Coords = double(tokens(:,1).');
      
    elseif isSyntax3
      % Initialize
      Coords = zeros(1,length(matched));
      Transitions = cell(1,nSpins);
      spinIndexPresent = false(1,length(matched));
      
      % String pattern to get tokens
      expr = '([exyz+-])\((\d+)\|?(\d+)?\)(\d?)';

      % loop over all comma separated entries
      for itoken = 1:length(matched)
        % special case: string does not request a transition/level 
        if ~any(matched{itoken}{1}=='(')
          % get tokens
          tokens = regexp(matched{itoken}{1},'([exyz+-])(\d?)','tokens');
          token = tokens{1};
          
          % get spin coordinate
          Coords(itoken) = token{1};
          % empty
          Transitions{itoken} = [];
          
          % Store spin index...
          if isempty(token{2})
            % ... from the counter...
            Spins(itoken) = itoken;
          else
            %... if a spin index was provided
            Spins(itoken) = str2double(token{2});
            spinIndexPresent(itoken) = true;
          end
        else
          % get tokens
          tokens = regexp(matched{itoken}{1},expr,'tokens');
          token = tokens{1};
          
          % get spin coordinate
          Coords(itoken) = token{1};
          
          % get spin index if available
          spinIndexPresent(itoken) = ~isempty(token{4});
          if ~spinIndexPresent(itoken)
            Spins(itoken) = itoken;
          else
            Spins(itoken) = str2double(token{4});
          end
          
          % get transition...
          if ~isempty(token{2})
            if ~isempty(token{3})
              Transitions{itoken} = cellfun(@str2num,token(2:3));
              % error if xyz+- is called connect two identical levels
              if ~isempty(regexp(token{1},'[xyz+-]','match')) && Transitions{itoken}(1) == Transitions{itoken}(2)
                message = ['sop: The component ''' matched{itoken}{1} ''' of your spin operator connects to identical levels.'];
                error(message)
              % error if e is called with two different components
              elseif ~isempty(regexp(token{1},'[e]','match')) && Transitions{itoken}(1) ~= Transitions{itoken}(2)
                message = ['sop: The component ''' matched{itoken}{1} ''' can not connect two different levels. Use e(L1) instead of e(L1|L2).'];
                error(message)
              end
          % ... or level
            else
              Transitions{itoken} = str2double(token{2});
              % error if xyz+- is called with only a level
              if ~isempty(regexp(token{1},'[xyz+-]','match'))
                message = ['sop: The component ''' matched{itoken}{1} ''' of your operator must connect two (different) levels.'];
                error(message)
              end
            end
          end
        end
      end
      
      % Further checking of Syntax
      % Are all components declared in the same syntax?
      if ~all(spinIndexPresent) && ~all(~spinIndexPresent)
        error('sop: Please use a consistent syntax to declare your spin operator ''%s''.' ,lower(varargin{1}))
      end
      % Is a component missing for a syntax of 'x,x'?
      if ~any(spinIndexPresent) && length(matched) ~= nSpins
        error('sop: Could not determine what ''%s'' is for the given spin system.',lower(varargin{1}));
      elseif all(~spinIndexPresent)
        Spins = 1:nSpins;
      end
      % Repeated spin
      if numel(unique(Spins))<numel(Spins)
        error('sop: Repeated spin in %s',lower(varargin{1})');
      end
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
  
  % Look for a transition selective operator
  if exist('Transitions','var') && ~isempty(Transitions{iSpin})
    TselectiveOp = true;
  else
    TselectiveOp = false;
  end
  
  % Component switchyard
  %----------------------------------------
  switch Comps(iSpin)
    case {'x',1} % x component
      if TselectiveOp
        ib = [Transitions{iSpin}(1); Transitions{iSpin}(2)];
        jb = [Transitions{iSpin}(2); Transitions{iSpin}(1)];
        sb = [0.5; 0.5];
      else
        m = (1:n-1).';
        ib = [m; m+1];
        jb = [m+1; m];
        Dia = 1/2*sqrt(m.*m(end:-1:1));
        sb = [Dia; Dia];
      end
    case {'y',2} % y component
      if TselectiveOp
        ib = [Transitions{iSpin}(1); Transitions{iSpin}(2)];
        jb = [Transitions{iSpin}(2); Transitions{iSpin}(1)];
        sb = [-0.5i; 0.5i];
      else
        m = (1:n-1).';
        Dia = -0.5i*sqrt(m.*m(end:-1:1));
        ib = [m; m+1];
        jb = [m+1; m];
        sb = [Dia; -Dia];
      end
    case {'z',3} % z component
      if TselectiveOp
        if length(Transitions{iSpin}) == 1
          ib = [Transitions{iSpin}(1)];
          jb = [Transitions{iSpin}(1)];
          sb = [1];
        else
          ib = [Transitions{iSpin}(1); Transitions{iSpin}(2)];
          jb = [Transitions{iSpin}(1); Transitions{iSpin}(2)];
          sb = [0.5; -0.5];
        end
      else
        m = (1:n).';
        ib = m;
        jb = m;
        sb = I+1-m;
      end
    case {'+',4} % up shift
      if TselectiveOp
        ib = [Transitions{iSpin}(1)];
        jb = [Transitions{iSpin}(2)];
        sb = [1];
      else
        m = (1:n-1).';
        ib = m;
        jb = m+1;
        sb = sqrt(m.*m(end:-1:1));
      end
    case {'-',5} % down shift
      if TselectiveOp
        ib = [Transitions{iSpin}(2)];
        jb = [Transitions{iSpin}(1)];
        sb = [1];
      else
        m = (1:n-1).';
        ib = m+1;
        jb = m;
        sb = sqrt(m.*m(end:-1:1));
      end
    case 'a' % alpha, for spin-1/2 only
      ib = 1;
      jb = 1;
      sb = 1;
    case 'b' % beta, for spin-1/2 only
      ib = 2;
      jb = 2;
      sb = 1;
    case {'e',0} % identity
      if TselectiveOp
        ib = [Transitions{iSpin}(1)];
        jb = [Transitions{iSpin}(1)];
        sb = [1];
      else
        m = (1:n).';
        ib = m;
        jb = m;
        sb = ones(n,1);
      end
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
