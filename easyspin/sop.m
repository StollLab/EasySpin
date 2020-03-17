% sop  Spin operator matrices
%
%   Op = sop(SpinSystem, Comps)
%   [Op1,Op2,...] = sop(SpinSystem,Comps1,Comps2,...)
%   ... = sop(...,'sparse')
%
%   Spin operator matrix of the spin system SpinSystem in the standard
%   |mS,mI,...> representation.
%
%   If more than one component is given, a matrix is computed for each component.
%
%   Input:
%    SpinSystem  vector of spin quantum numbers
%                 or a spin system specification structure
%    Comps       (a) string specifying the operator, in several possible ways:                 
%                    - specify component 'e','x','y','z','+','-' for each spin,
%                      indicating E,Sx,Sy,Sz,S+,S-
%                    - specify component and spin index, e.g. 'x2,z3'
%                    - specify transition after component, e.g. 'x(1|3)' or
%                      'x(1|3)2,z3' or '+(1|2)1,e(3)2'
%                (b) numeric array, with each row giving [i c], with spin index
%                    i and component index c (c=1 is 'x', 2 is 'y', 3 is 'z',
%                    4 is '+', 5 is '-', 0 is 'e');
%
%   Output:
%    Op          spin operator matrix as requested
%
%   Examples:
%     SxIy = sop([1/2 1],'xy')    % returns SxIy for a S=1/2, I=1 system
%
%     SeIp = sop([1/2 1/2],'e+')  % returns SeI+ for a S=I=1/2 system
%
%     [Sx,Sy,Sz] = sop(1/2,'x','y','z')  % computes three matrices in one call
%
%     Sxc = sop(5/2,'x(3|4)') % Sx on central transition -1/2<->+1/2

function varargout = sop(SpinSystem,varargin)

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

if nargin==0, help(mfilename); return; end

% Check first input
if isstruct(SpinSystem)
  SpinVec = spinvec(SpinSystem);
else
  SpinVec = SpinSystem;
end
if ~isnumeric(SpinVec) || ~isvector(SpinVec)
  error('First input must be a 1D array of spin quantum numbers.');
end
if any(SpinVec<1/2) || any(mod(SpinVec,1/2)) || ~isreal(SpinVec)
  error('First input must contain valid spin quantum numbers (1/2, 1, 3/2, etc).');
end
nSpins = numel(SpinVec);

if nargin==1
  error('Not enough input arguments - at least 2 are needed!');
end

sparseOutput = ischar(varargin{end}) && strcmpi(varargin{end},'sparse');

if sparseOutput
  OperatorSpec = varargin(1:end-1);
else
  OperatorSpec = varargin;
end

if isnumeric(OperatorSpec{1})
  
  % Operator specification is a numeric array
  
  if numel(OperatorSpec)~=1
    error('Incorrect number of input arguments.')
  end
  a = OperatorSpec{1};
  if size(a,2)~=2
    error('For numeric input, the array specifying spins and components must be Nx2.');
  end
  Spins = a(:,1);
  Coords = a(:,2);
  
else
  
  % Operator specification is a character array
  
  if nargout~=numel(OperatorSpec) && nargout>0
    error('sop: Number of output arguments (%d) and number of requested operator matrices (%d) do not match.',nargout,numel(OperatorSpec));
  end
  
  if numel(OperatorSpec)>1
    
    for k = 1:numel(OperatorSpec)
      if sparseOutput
        varargout{k} = sop(SpinVec,OperatorSpec{k},'sparse');
      else
        varargout{k} = sop(SpinVec,OperatorSpec{k});
      end
    end
    return
    
  else
    
    OperatorSpec = varargin{1};
    
    % Determine which type of syntax is used in the operator specification
    %---------------------------------------------------------------------------
    syntax1regexp = sprintf('^[exyzab+-]{%d}$',nSpins);
    isSyntax1 = ~isempty(regexp(OperatorSpec,syntax1regexp,'match'));
    if ~isSyntax1
      syntax2regexp = sprintf('^([exyzab+-]\\d+){1,%d}$',nSpins);
      isSyntax2 = ~isempty(regexp(OperatorSpec,syntax2regexp,'match'));
      if ~isSyntax2
        OperatorSpec = regexp(OperatorSpec,',','split');
        syntax3regexp = '^[exyzab+-](\(\d+(\|\d+)?\))?(\d+)?$';
        matched = regexp(OperatorSpec,syntax3regexp,'match');
        isSyntax3 = ~any(cellfun(@isempty, matched));
        if ~isSyntax3
          error('sop: Cannot parse ''%s'' for the given spin system.',varargin{1});
        end
      end
    end
    
    % Parse inputs depending on type of syntax
    %---------------------------------------------------------------------------
    if isSyntax1
      
      Coords = OperatorSpec;
      Spins = 1:nSpins;
      
    elseif isSyntax2
      
      % Get tokens from pattern
      pattern = '([exyzab+-]\d+)+?';
      tokens = char(regexp(OperatorSpec,pattern,'match'));
      % Get list of spin indices
      Spins = str2num(tokens(:,2:end)).'; %#ok<ST2NM>
      if any(Spins>nSpins)
        error('sop: The spin system only contains %d spins, but you requested ''%s''.',nSpins,varargin{1});
      end
      if numel(unique(Spins))<numel(Spins)
        error('sop: Duplicate spin index in ''%s''.',OperatorSpec);
      end
      % Get list of requested components
      Coords = tokens(:,1).';
      
    elseif isSyntax3
      
      % Initialize
      Coords = '';
      Transitions = cell(1,nSpins);
      spinIndexPresent = false(1,length(matched));
      
      % Regexp to get tokens
      expr = '([exyzab+-])\((\d+)\|?(\d+)?\)(\d?)';

      % Loop over all comma separated entries
      for itoken = 1:length(matched)
        
        str = matched{itoken}{1};
        
        if ~any(str=='(')
          % Case 1: string does not request a transition/level 
          
          % Get tokens
          tokens = regexp(str,'([exyzab+-])(\d?)','tokens');
          token = tokens{1};
          
          % Get operator components
          Coords(itoken) = token{1};
          % No transition information
          Transitions{itoken} = [];
          
          % Store spin index...
          spinIndexPresent(itoken) = ~isempty(token{2});
          if ~spinIndexPresent(itoken)
            % ... from the counter.
            Spins(itoken) = itoken;
          else
            %... if a spin index was provided.
            Spins(itoken) = str2double(token{2});
          end
          
        else
          % Case 2: String contains transition/level syntax (|)
          
          % Get tokens
          tokens = regexp(str,expr,'tokens');
          token = tokens{1};
          
          % Get operator component
          Coords(itoken) = token{1};
          if any(token{1}=='ab')
            error('sop: ''a'' or ''b'' cannot be used in conjunction with levels in ''%s''',varargin{1});
          end
          
          % Get spin index if available
          spinIndexPresent(itoken) = ~isempty(token{4});
          if ~spinIndexPresent(itoken)
            Spins(itoken) = itoken;
          else
            Spins(itoken) = str2double(token{4});
          end
          
          % Get transition/level information
          if ~isempty(token{2})
            if ~isempty(token{3})
              Transitions{itoken} = cellfun(@str2num,token(2:3));
              % Error if xyz+- is connecting two identical levels
              if ~isempty(regexp(token{1},'[xyz+-]','match')) && Transitions{itoken}(1) == Transitions{itoken}(2)
                message = ['sop: The component ''' str ''' of your spin operator connects identical levels.'];
                error(message);
              % Error if e is called with two different components
              elseif ~isempty(regexp(token{1},'[e]','match')) && Transitions{itoken}(1) ~= Transitions{itoken}(2)
                message = ['sop: The component ''' str ''' can not connect two different levels. Use e(L1) instead of e(L1|L2).'];
                error(message);
              % Error if the two level indices are not in ascending order
              elseif diff(Transitions{itoken})<0
                message = ['sop: The two levels (L1|L2) in ''' str ''' must be ascending order such that L1<=L2.'];
                error(message);
              end
            else
              Transitions{itoken} = str2double(token{2});
              % error if xyz+- is called with only a level
              if ~isempty(regexp(token{1},'[xyz+-]','match'))
                message = ['sop: The component ''' str ''' of your operator must connect two (different) levels.'];
                error(message);
              end
            end
          end
          
        end
      end
      
      % Further syntax checking
      % Make sure either all or no components have explicit spin indices
      if ~all(spinIndexPresent) && ~all(~spinIndexPresent)
        error('sop: Please give spin indices for all components, or omit them all. ''%s'' is inconsistent.' ,varargin{1})
      end
      % Is a component missing for a syntax of 'x,x'?
      if ~any(spinIndexPresent) && length(matched) ~= nSpins
        error('sop: Could not determine what ''%s'' is for the given spin system.',varargin{1});
      elseif all(~spinIndexPresent)
        Spins = 1:nSpins;
      end
      % Assert that no spin index occurs more than once
      if numel(unique(Spins))<numel(Spins)
        error('sop: Repeated spin index in ''%s''.',varargin{1});
      end
    end
    
    if any(SpinVec(Coords=='a')~=1/2) || ...
       any(SpinVec(Coords=='b')~=1/2)
      error('''a'' and ''b'' work only for spin-1/2.');
    end
    
  end
end

% Initialize with identity operator for each spin
Components = repmat('e',1,nSpins);
% Assign specified components for specified spins
Components(Spins) = Coords;

% Set starting operator matrix
Op = sparse(1);

% Run over all spins
for iSpin = 1:nSpins
  I = SpinVec(iSpin);
  if I<=0, continue; end
  
  n = 2*I+1;
  
  % Look for a transition selective operator
  TselectiveOp = exist('Transitions','var') && ~isempty(Transitions{iSpin});
  
  % Component switchyard
  %----------------------------------------
  switch Components(iSpin)
    case {'x',1} % x component
      if TselectiveOp
        r = [Transitions{iSpin}(1); Transitions{iSpin}(2)];
        c = [Transitions{iSpin}(2); Transitions{iSpin}(1)];
        val = [0.5; 0.5];
      else
        m = 1:n-1;
        r = [m; m+1];
        c = [m+1; m];
        dia = 1/2*sqrt(m.*m(end:-1:1));
        val = [dia; dia];
      end
    case {'y',2} % y component
      if TselectiveOp
        r = [Transitions{iSpin}(1); Transitions{iSpin}(2)];
        c = [Transitions{iSpin}(2); Transitions{iSpin}(1)];
        val = [-0.5i; 0.5i];
      else
        m = 1:n-1;
        dia = -0.5i*sqrt(m.*m(end:-1:1));
        r = [m; m+1];
        c = [m+1; m];
        val = [dia; -dia];
      end
    case {'z',3} % z component
      if TselectiveOp
        r = [Transitions{iSpin}(1); Transitions{iSpin}(2)];
        c = [Transitions{iSpin}(1); Transitions{iSpin}(2)];
        val = [0.5; -0.5];
      else
        m = 1:n;
        r = m;
        c = m;
        val = I+1-m;
      end
    case {'+',4} % up shift
      if TselectiveOp
        r = Transitions{iSpin}(1);
        c = Transitions{iSpin}(2);
        val = 1;
      else
        m = 1:n-1;
        r = m;
        c = m+1;
        val = sqrt(m.*m(end:-1:1));
      end
    case {'-',5} % down shift
      if TselectiveOp
        r = Transitions{iSpin}(2);
        c = Transitions{iSpin}(1);
        val = 1;
      else
        m = 1:n-1;
        r = m+1;
        c = m;
        val = sqrt(m.*m(end:-1:1));
      end
    case 'a' % alpha, for spin-1/2 only
      r = 1;
      c = 1;
      val = 1;
    case 'b' % beta, for spin-1/2 only
      r = 2;
      c = 2;
      val = 1;
    case {'e',0} % identity
      if TselectiveOp
        r = Transitions{iSpin}(1);
        c = Transitions{iSpin}(1);
        val = 1;
      else
        m = 1:n;
        r = m;
        c = m;
        val = ones(n,1);
      end
    otherwise
      error('Unknown operator specification.');
  end
  
  M_ = sparse(r,c,val,n,n);
  Op = kron(Op,M_);
  
end

% Convert sparse to full matrix if required
if ~sparseOutput
  Op = full(Op);
end

varargout = {Op};

return
