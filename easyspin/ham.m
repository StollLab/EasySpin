% ham  Spin Hamiltonian 
%
%   [F,GxM,GyM,GzM] = ham(Sys)
%   [F,GzL] = ham(Sys,B0)
%   H = ham(Sys,B0)
%   H = ham(Sys)
%   __ = ham(__,'sparse')
%
%   Constructs the spin Hamiltonian matrix, or its field-independent and
%   field-dependent components, for the spin system Sys.
%
%   Input:
%     Sys       spin system specification structure
%     B0        vector specifying the static magnetic field (in mT) in the
%               molecular frame (e.g. [350 0 0] is along the molecular x axis)
%     'sparse'  if given, sparse instead of full matrices are returned
%
%   Output:
%     H            complete spin Hamiltonian (MHz)
%     F            field-independent part of spin Hamiltonian (MHz)
%     GxM,GyM,GzM  field-dependent spin Hamiltonian components (MHz/mT)
%                  x, y, z axes of molecular frame
%     GzL          field-dependent spin Hamiltonian (MHz/mT), with field
%                  along z axis of lab frame
%
%   The spin Hamiltonian components are defined such that
%           H = F + B0(1)*GxM + B0(2)*GyM + B0(3)*GzM

function varargout = ham(SpinSystem,B0,opt)

if nargin==0, help(mfilename); return; end

if nargin<2, B0 = []; end
if nargin<3, opt = ''; end
if ~ischar(opt)
  error('Third argument must be a string, ''sparse'' or ''''.');
end
sparseMatrices = strcmp(opt,'sparse');
if sparseMatrices
  sp = 'sparse';
else
  sp = 'full';
end

[Sys,err] = validatespinsys(SpinSystem);
error(err);

if ~isempty(B0)
  if ~isnumeric(B0) || numel(B0)~=3 || ~isreal(B0)
    error('Magnetic field vector must be 3-element array!');
  end
end

% Build field-independent part of Hamiltonian
F = ham_zf(Sys,[],sp) + ham_ee(Sys,[],sp) + ...
    ham_hf(Sys,[],[],sp) + ham_nq(Sys,[],sp) + ham_nn(Sys,[],sp) + ...
    ham_so(Sys,[],sp) + ham_cf(Sys,[],sp);

% Build field-dependent part of Hamiltonian
[GxM,GyM,GzM] = ham_ez(Sys,[],sp);
if Sys.nNuclei>0
  [GxMn,GyMn,GzMn] = ham_nz(Sys,[],sp);
  GxM = GxM + GxMn;
  GyM = GyM + GyMn;
  GzM = GzM + GzMn;
end
if Sys.nL>0
  [GxMo,GyMo,GzMo] = ham_oz(Sys,[],sp);
  GxM = GxM + GxMo;
  GyM = GyM + GyMo;
  GzM = GzM + GzMo;
end

% Look for higher-order Zeeman terms
sysfields = fieldnames(Sys);
highest = 0;
higherOrder = false;
higherzeeman = strncmp(sysfields,'Ham',3).';
if any(higherzeeman) 
  for n = find(higherzeeman)
    if any(Sys.(sysfields{n}))
      higherOrder = true;
      order = str2double(sysfields{n}(4));
      if order > highest && order < 4, highest = order; end
    end
  end
end

if isempty(B0)
  
  if higherOrder
    % full tensors up to the highest used order will be provided
    zHo = ham_ezho(Sys,[],opt);
    switch highest
      case 0
        if nargout==4
          varargout = {F+zHo,GxM,GyM,GzM};
        elseif nargout==1
          varargout{1} = {F+zHo,GxM,GyM,GzM};
        else
          error('1 or 4 output arguments expected!');
        end
      case 1
        if nargout==4
          varargout = {F+zHo{1},GxM+zHo{2}{1},GyM+zHo{2}{2},GzM+zHo{2}{3}};
        elseif nargout==1
          varargout{1} = {F+zHo{1},GxM+zHo{2}{1},GyM+zHo{2}{2},GzM+zHo{2}{3}};
        else
          error('1 or 4 output arguments expected!');
        end
      otherwise
        zHo{1} = zHo{1}+F;
        Gn = {GxM,GyM,GzM};
        for k = 3:-1:1
          zHo{2}{k} = zHo{2}{k}+Gn{k};
        end
        if nargout==highest+1
          varargout = zHo;
        elseif nargout==1
          varargout{1} = zHo;
        else
          error('Wrong number of output arguments!');
        end
    end
  else
    
    if nargout==4
      varargout = {F,GxM,GyM,GzM};
    elseif nargout==1
      varargout{1} = {{F,GxM,GyM,GzM}};
    else
      error('1 or 4 output arguments expected!');
    end
  end
  
  return
  
end


% Field is provided
%------------------------------------------------------------------------
if norm(B0)>0
  nB0 = B0/norm(B0);
else
  nB0 = [0 0 0];
end
GzL = nB0(1)*GxM + nB0(2)*GyM + nB0(3)*GzM;

if higherOrder
  % ~isempty(B0) && higherOrder
  if nargout==1 || nargout==0
    varargout = {F + norm(B0)*GzL + ham_ezho(Sys,B0,[],opt)};
  elseif nargout==2 && highest<2
    zHo = ham_ezho(Sys,[],opt);
    if highest == 1
      GzL = GzL + nB0(1)*zHo{2}{1} + nB0(2)*zHo{2}{2} + nB0(3)*zHo{2}{3};
      F = F + zHo{1};
    else
      F = F +zHo;
    end
    varargout = {F,GzL};
  else
    error('Wrong number of output arguments!');
  end
else
  % ~isempty(B0) && ~higherOrder
  if nargout==1 || nargout==0
    varargout = {F + norm(B0)*GzL};
  elseif nargout==2
    varargout = {F,GzL};
  else
    error('Wrong number of output arguments!');
  end
end

end
