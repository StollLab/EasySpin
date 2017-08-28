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

if OldSyntax
  
  Spins = varargin{1};
  Coords = varargin{2};
  
else
  
  if nargout~=numel(Coords)
    error('Number of output arguments (%d) and number of requested operator matrices (%d) do not match.',nargout,numel(Coords));
  end
  
  if numel(Coords)>1
    for k = 1:numel(Coords)
      if SparseOutput
        varargout{k} = sop(SpinVec,Coords{k},'sparse');
      else
        varargout{k} = sop(SpinVec,Coords{k});
      end
    end
    return
  else
    
    Coords = lower(varargin{1});
    
    syntax1 = sprintf('^[exyz+\\-pmab]{%d}$',nSpins);
    isSyntax1 = ~isempty(regexp(Coords,syntax1,'match'));
    if ~isSyntax1
      syntax2 = sprintf('^([exyz+\\-pmab]\\d+){1,%d}$',nSpins);
      isSyntax2 = ~isempty(regexp(Coords,syntax2,'match'));
      if ~isSyntax2
        Coords = regexp(Coords,',','split');
        syntax3 = '^[exyz+\-pmab](\(\d+(\|\d+)?\))?(\d+)?$';
        matched = (regexp(Coords,syntax3,'match'));
        isSyntax3 = length(find(~cellfun('isempty', matched))) == length(Coords);
        if ~isSyntax3
          error('Could not determine what ''%s'' is for the given spin system.',lower(varargin{1}));
        end
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
    elseif isSyntax3
      % Get tokens and set up data structure
      Syntaxtypes = zeros(1,length(matched));
      Coords = Syntaxtypes;
      Spins = Syntaxtypes;
      Transitions = cell(1,nSpins);
      % Loop over all the comma separated operators, e.g.
      % x(1|2)1,x(1|2)2/x(1|2),x(1|2) or x1,x2/x,x
      for itoken = 1 : length(matched)
       token = char(matched{itoken}); 
       % split at paranthesis if possible - for transitions selection
       splitted = regexp(token,'[()]','split');
       Coords(itoken) = double(splitted{1}(1));
       % e.g. length(splitted) > 1: x(1|2)1 or x(1|2)
       if length(splitted) > 1
         if isempty(splitted{end})
           % x(1|2)
           Syntaxtypes(itoken) = 1;
           % Get spin index temporarily for indexing for the selective 
           % transtion 
           Spin = itoken;
         else
           % x(1|2)
           Syntaxtypes(itoken) = 2;
           % Store spin number in vector for processing later ...
           Spins(itoken) = str2double(splitted{end});
           % ... and temporarily for indexing for the selective transtion  
           Spin = Spins(itoken);
         end
         % Find requested level(s)
         Selector = regexp(splitted{2},'\|','split');
         % Populationselective, e.g  by e(1)1/e(1) or z(1)/z(1)1
         if length(Selector) == 1
           if isempty(regexp(splitted{1}(1),'[ez]','match'))
             error('sop: You need to provide the second connected level for your operator or change the operator type.')
           end
           Transitions{Spin} = str2double(Selector);
         % Populationselective, e.g  e(1|1)1/e(1|1) or z(1|1)/z(1|1)1
         elseif Selector{1} == Selector{2}
           if ~isempty(regexp(splitted{1}(1),'[xyabpm+-]','match'))
             error('sop: The transition selective operator you selected connects to identical levels.')
           end
            Transitions{Spin} = str2double(Selector{1});
         % Operators that connect two leves, e.g x(1|2)1/x(1|2)
         else
           if ~isempty(regexp(splitted{1}(1),'[e]','match'))
             error('sop: Only one level can be selected for a population.')
           end
           Transitions{Spin}(1) = str2double(Selector{1});
           Transitions{Spin}(2) = str2double(Selector{2});
         end
       else
         % length(splitted) = 1: x1 or x
         if length(splitted{1}) == 1
           Syntaxtypes(itoken) = 1;
         else
           Syntaxtypes(itoken) = 2;
           Spins(itoken) = str2double(splitted{1}(2:end));
         end
       end
      end
      % Further checking of Syntax
      if numel(unique(Syntaxtypes)) > 1
        error('sop: Please use consistent syntax to declare your spin operator.')
      end
      if ~any(Syntaxtypes == 2) && length(matched) ~= nSpins
        error('Could not determine what ''%s'' is for the given spin system.',lower(varargin{1}));
      elseif unique(Syntaxtypes) == 1
        Spins = 1:nSpins;
      end
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
    case {'+','p',4} % up shift
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
    case {'-','m',5} % down shift
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
