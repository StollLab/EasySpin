% sham  Spin Hamiltonian 
%
%   [F,Gx,Gy,Gz] = sham(sys)
%   H = sham(sys, B)
%   ... = sham(sys, B,'sparse')
%
%   Constructs a spin Hamiltonian.
%
%   Input:
%   - 'sys': spin system specification structure
%   - 'B': 1x3 vector specifying the magnetic field [mT]
%   - 'sparse': if present, sparse instead of full matrices are returned
%
%   Output:
%   - 'H': complete spin Hamiltonian [MHz]
%   - 'F',Gx','Gy','Gz' ([MHz] and [MHz/mT])
%      H = F + B(1)*Gx + B(2)*Gy + B(3)*Gz

function varargout = sham(SpinSystem,B0,opt)

if (nargin==0), help(mfilename); return; end

if nargin<2, B0 = []; end
if nargin<3, opt = ''; end
if ~ischar(opt)
  error('Third argument must be a string, ''sparse''.');
end
sparseResult = strcmp(opt,'sparse');

[Sys,err] = validatespinsys(SpinSystem);
error(err);


sysfields = fieldnames(Sys);
highest = 0;
higherOrder = false;
higherzeeman = strncmp(sysfields,'ZB',2).';
if any(higherzeeman) 
  for n=find(higherzeeman)
    if isfield(Sys.(sysfields{n}), 'vals') 
      if ~iscell(Sys.(sysfields{n}).vals) && any(Sys.(sysfields{n}).vals(:))
        higherOrder = true;
        order = str2num(sysfields{n}(3));
        if order > highest && order < 4, highest = order;end
      else
        t = 0;
        for m= length(Sys.(sysfields{n}).vals)
          t = t + any(Sys.(sysfields{n}).vals{m}(:));
        end
        if t
          higherOrder = true;
          order = str2num(sysfields{n}(3));
          if order > highest && order < 4, highest = order;end
        end
      end
    end
  end
end


% Field-independent interactions: ZFI, NQI, HFI, EEI
% Zeeman interaction: EZI, NZI
if sparseResult
  F = zfield(Sys,[],'sparse') + nquad(Sys,[],'sparse') + ...
    hfine(Sys,[],'sparse') + eeint(Sys,[],'sparse');
  [GxM,GyM,GzM] = zeeman(Sys,[],'sparse');
else
  F = zfield(Sys) + nquad(Sys) + hfine(Sys) + eeint(Sys);
  [GxM,GyM,GzM] = zeeman(Sys);
end


if higherOrder
  if isempty(B0)
    %full tensors up to the highest used orer will be provided
    zHo = zeemanho(Sys,[],opt);
    switch highest
      case 0
        if (nargout==4)
          varargout = {F+zHo,GxM,GyM,GzM};
        elseif nargout==1
          varargout{1} = {F+zHo,GxM,GyM,GzM};
        else
          error('1 or 4 output arguments expected!');
        end
      case 1
        if (nargout==4)
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
        if(nargout==highest+1)
          varargout =zHo;
        elseif nargout==1
          varargout{1}=zHo;
        else
          error('Wrong number of output arguments!');
        end
    end
  else
    if numel(B0)~=3
      error('Magnetic field must be 3-element vector!');
    end
    if norm(B0)>0
      nB0 = B0/norm(B0);
    else
      nB0 = [0 0 0];
    end
    GzL = nB0(1)*GxM + nB0(2)*GyM + nB0(3)*GzM;    
    if (nargout==1) || (nargout==0)
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
  end
else
  % arrange the output
  if isempty(B0)
    
    if (nargout==4)
      varargout = {F,GxM,GyM,GzM};
    elseif nargout==1
      varargout{1} = {F,GxM,GyM,GzM};
    else
      error('1 or 4 output arguments expected!');
    end
    
  else
    
    if numel(B0)~=3
      error('Magnetic field must be 3-element vector!');
    end
    if norm(B0)>0
      nB0 = B0/norm(B0);
    else
      nB0 = [0 0 0];
    end
    GzL = nB0(1)*GxM + nB0(2)*GyM + nB0(3)*GzM;
    if (nargout==1) || (nargout==0)
      varargout = {F + norm(B0)*GzL};
    elseif nargout==2
      varargout = {F,GzL};
    else
      error('Wrong number of output arguments!');
    end
    
  end
end
return
