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

% Field-independent interactions: ZFI, NQI, HFI, EEI
% Zeeman interaction: EZI, NZI
if sparseResult
  F = zfield(Sys,[],'sparse') + nquad(Sys,[],'sparse') + ...
    hfine(Sys,[],'sparse') + eeint(Sys,'sparse');
  [GxM,GyM,GzM] = zeeman(Sys,[],'sparse');
else
  F = zfield(Sys) + nquad(Sys) + hfine(Sys) + eeint(Sys);
  [GxM,GyM,GzM] = zeeman(Sys);
end


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

return
