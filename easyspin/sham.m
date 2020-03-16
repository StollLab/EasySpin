% sham  Spin Hamiltonian 
%
%   [F,GxM,GyM,GzM] = sham(Sys)
%   [F,GzL] = sham(Sys,B0)
%   H = sham(Sys,B0)
%   H = sham(Sys)
%   __ = sham(__,'sparse')
%
%   Constructs the spin Hamiltonian, or its field-independent and field-
%   dependent components, for the spin system Sys.
%
%   Input:
%     sys       spin system specification structure
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
%           GzL = (B0(1)*GxM+B0(2)*GyM+B0(3)*GzM)/norm(B0)
%           H = F + norm(B0)*GzL

function varargout = sham(SpinSystem,B0,opt)

if nargin==0, help(mfilename); return; end

if nargin<2, B0 = []; end
if nargin<3, opt = ''; end
if ~ischar(opt)
  error('Third argument must be a string, ''sparse''.');
end
sparseMatrices = strcmp(opt,'sparse');

[Sys,err] = validatespinsys(SpinSystem);
error(err);

if ~isempty(B0)
  if ~isnumeric(B0) || numel(B0)~=3 || ~isreal(B0)
    error('Magnetic field vector must be 3-element array!');
  end
end


% Field-independent interactions: ZFI, NQI, HFI, EEI, NNI
% Zeeman interaction: EZI, NZI
if sparseMatrices
  sp = 'sparse';
else
  sp = 'full';
end

F = zfield(Sys,[],sp) + eeint(Sys,[],sp) + ...
    hfine(Sys,[],sp)  + nquad(Sys,[],sp) + nnint(Sys,[],sp) + ...
    soint(Sys,[],sp)  + crystalfield(Sys,[],sp);
[GxM,GyM,GzM] = zeeman(Sys,[],sp);


sysfields = fieldnames(Sys);
highest = 0;
higherOrder = false;
higherzeeman = strncmp(sysfields,'Ham',3).';
if any(higherzeeman) 
  for n = find(higherzeeman)
    if any(Sys.(sysfields{n}))
      higherOrder = true;
      order = str2double(sysfields{n}(4));
      if order > highest && order < 4, highest = order;end
    end
  end
end

if isempty(B0)
  
  if higherOrder
    % full tensors up to the highest used orer will be provided
    zHo = zeemanho(Sys,[],opt);
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


if higherOrder
  % ~isempty(B0) && higherOrder
  if norm(B0)>0
    nB0 = B0/norm(B0);
  else
    nB0 = [0 0 0];
  end
  GzL = nB0(1)*GxM + nB0(2)*GyM + nB0(3)*GzM;
  if nargout==1 || nargout==0
    varargout = {F + norm(B0)*GzL + zeemanho(Sys,B0,[],opt)};
  elseif nargout==2 && highest<2
    zHo = zeemanho(Sys,[],opt);
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
  if norm(B0)>0
    nB0 = B0/norm(B0);
  else
    nB0 = [0 0 0];
  end
  GzL = nB0(1)*GxM + nB0(2)*GyM + nB0(3)*GzM;
  if nargout==1 || nargout==0
    varargout = {F + norm(B0)*GzL};
  elseif nargout==2
    varargout = {F,GzL};
  else
    error('Wrong number of output arguments!');
  end
end

return
