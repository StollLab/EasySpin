% sham  Spin Hamiltonian 
%
%   [F,Gx,Gy,Gz] = sham(sys)
%   [F,G] = sham(sys, B)
%   H = sham(sys, B)
%
%   Constructs a spin Hamiltonian.
%
%   Input:
%   - 'sys': spin system specification structure
%   - 'B': 1x3 vector specifying the magnetic field [mT]
%
%   Output:
%   - 'H': complete spin Hamiltonian [MHz]
%   - 'F','G',Gx','Gy','Gz' ([MHz] and [MHz/mT])
%      H = F + B(1)*Gx + B(2)*Gy + B(3)*Gz = F + |B|*G

function varargout = sham(SpinSystem,B0)

if (nargin==0), help(mfilename); return; end

[Sys,err] = validatespinsys(SpinSystem);
error(err);

% Field-independent interactions: ZFI, NQI, HFI, EEI
F = zfield(Sys) + nquad(Sys) + hfine(Sys) + eeint(Sys);

% Zeeman interaction: EZI, NZI
[GxM,GyM,GzM] = zeeman(Sys);

% arrange the output
switch nargin
case 1
  if (nargout==4)
    varargout = {F,GxM,GyM,GzM};
  elseif nargout==1
    varargout{1} = {F,GxM,GyM,GzM};
  else
    error('4 output arguments expected!');
  end
  
case 2
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
