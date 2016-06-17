function [err,data] = test(opt,olddata)

% Check excitation profiles calculated by pulse()
%--------------------------------------------------------------------------

% Rectangular pulses
Exp.tp = [0.016 0.032]; % us
Exp.Flip = [pi/2 pi];
Exp.Phase = [pi/2 pi/2];

Opt.Detect = 'all';
[~,~,p] = pulse(Exp,1:2,Opt);

% Calculate inversion profile
Sx = sop(1/2,'x');
Sy = sop(1/2,'y');
Sz = sop(1/2,'z');

Mx(1:2,1:numel(p.offsets)) = 0;
My(1:2,1:numel(p.offsets)) = 0;
Mz(1:2,1:numel(p.offsets)) = 0;
for k = 1:2
  
  Amplitude = (Exp.Flip(k)/Exp.tp(k))/(2*pi);
  rho0 = -Sz;
  
  for i = 1:numel(p.offsets)
    
    H = p(k).offsets(i)*Sz + Amplitude*Sy;
    U = expm(-2i*pi*H*Exp.tp(k));
    rho = U*rho0*U';
    
    Mx(k,i) = -2*trace(Sx*rho);
    My(k,i) = -2*trace(Sy*rho);
    Mz(k,i) = -2*trace(Sz*rho);
    
  end
  
end

suberr(1) = ~areequal(Mx(1,:),p(1).Mx,1e-12);
suberr(2) = ~areequal(Mx(2,:),p(2).Mx,1e-12);
suberr(3) = ~areequal(My(1,:),p(1).My,1e-12);
suberr(4) = ~areequal(My(2,:),p(2).My,1e-12);
suberr(5) = ~areequal(Mz(1,:),p(1).Mz,1e-12);
suberr(6) = ~areequal(Mz(2,:),p(2).Mz,1e-12);

err(1) = any(suberr);

% 1st order sech/tanh with frequency offset
clearvars -except err
Exp.tp = 0.200; % us
Exp.PulseShape.Type = 'sech/tanh';
Exp.PulseShape.BW = 100; % MHz
Exp.PulseShape.beta = 15;
Exp.PulseShape.SweepDirection = -1;
Exp.Flip = pi;
Exp.TimeStep = 0.0005; % us
Exp.CenterFreq = 60; % MHz

t0 = 0:Exp.TimeStep:Exp.tp;
% Calculate pulse amplitude from flip angle
Qcrit = 8;
BW = Exp.PulseShape.BW/tanh(Exp.PulseShape.beta/2);
Amplitude = sqrt((Exp.PulseShape.beta*BW*Qcrit)/(2*pi));
% Amplitude modulation: sech
A = sech((Exp.PulseShape.beta/Exp.tp)*(t0-Exp.tp/2));
% Phase modulation
phi = (Exp.PulseShape.BW/(2*tanh(Exp.PulseShape.beta/2)))*(Exp.tp/Exp.PulseShape.beta)*log(cosh((Exp.PulseShape.beta/Exp.tp)*(t0-Exp.tp/2)));
phi = 2*pi*Exp.PulseShape.SweepDirection*phi;
% Pulse
y0 = Amplitude*A.*exp(1i*(phi+2*pi*Exp.CenterFreq*t0));

[~,~,p] = pulse(Exp);

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

err(2) = ~areequal(Mz,p.Mz,1e-3);

err = any(err);

data = [];