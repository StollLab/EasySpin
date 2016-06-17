function [err,data] = test(opt,olddata)

% Check excitation profiles calculated by pulse() with the BCH
% approximation
%--------------------------------------------------------------------------

% Sinc pulse
Exp.tp = 0.200; % us
Exp.Flip = pi;
Exp.PulseShape.Type = 'sinc';
Exp.PulseShape.zerocross = 0.050;

Opt.Detect = 'Sz';
Opt.Offsets = -100:1:100;
Opt.nBCH = 3;
[t,y,p] = pulse(Exp,Opt);

% Calculate inversion profile
Sx = sop(1/2,'x');
Sy = sop(1/2,'y');
Sz = sop(1/2,'z');

Mz(1:numel(p.offsets)) = 0;

rho0 = -Sz;

for i = 1:numel(p.offsets)
  
  rho = rho0;
  
  for j = 1:numel(t)
    H = p.offsets(i)*Sz + real(y(j))*Sx + imag(y(j))*Sy;
    U = expm(-2i*pi*H*(t(2)-t(1)));
    rho = U*rho*U';
  end
  
  Mz(i) = -2*trace(Sz*rho);
  
end
  
err(1) = ~areequal(Mz,p.Mz,0.5e-1);

% 1st order sech/tanh
clearvars -except err
Exp.tp = 0.100; % us
Exp.PulseShape.Type = 'sech/tanh';
Exp.PulseShape.BW = 60; % MHz
Exp.PulseShape.beta = 8;
Exp.Flip = pi;
Exp.TimeStep = 0.0005; % us

t0 = 0:Exp.TimeStep:Exp.tp;
% Calculate pulse amplitude from flip angle
Qcrit = 8;
BW = Exp.PulseShape.BW/tanh(Exp.PulseShape.beta/2);
Amplitude = sqrt((Exp.PulseShape.beta*BW*Qcrit)/(2*pi));
% Amplitude modulation: sech
A = sech((Exp.PulseShape.beta/Exp.tp)*(t0-Exp.tp/2));
% Phase modulation
phi = (Exp.PulseShape.BW/(2*tanh(Exp.PulseShape.beta/2)))*(Exp.tp/Exp.PulseShape.beta)*log(cosh((Exp.PulseShape.beta/Exp.tp)*(t0-Exp.tp/2)));
% Pulse
y0 = Amplitude*A.*exp(2i*pi*phi);

Opt.nBCH = 3;
Opt.Offsets = -120:1:120;
[t,y,p] = pulse(Exp,Opt);

% Inversion profile
Sx = sop(1/2,'x');
Sy = sop(1/2,'y');
Sz = sop(1/2,'z');

rho0 = -Sz;
Mz(1:numel(p.offsets)) = 0;
for i = 1:numel(p.offsets)
  
  H0 = p.offsets(i)*Sz;
  rho = rho0;
  
  for j = 1:numel(t0)
    H = H0 + real(y0(j))*Sx + imag(y0(j))*Sy;
    U = expm(-2i*pi*H*Exp.TimeStep);
    rho = U*rho*U';
  end 
  
  Mz(i) = -2*trace(Sz*rho);
  
end

err(2) = ~areequal(Mz,p.Mz,1e-2);

err = any(err);

data = [];