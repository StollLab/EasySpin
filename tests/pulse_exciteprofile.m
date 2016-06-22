function [err,data] = test(opt,olddata)

% Check excitation profiles calculated by pulse()
%--------------------------------------------------------------------------

% Rectangular pulses
Params(1).Flip = pi/2;
Params(1).tp = 0.016; % us
Params(1).Phase = pi/2;
Params(2).Flip = pi;
Params(2).tp = 0.032; % us
Params(2).Phase = pi/2;

Opt.Detect = 'all';
[~,~,p(1)] = pulse(Params(1),Opt);
[~,~,p(2)] = pulse(Params(2),Opt);

% Calculate inversion profile
Sx = sop(1/2,'x');
Sy = sop(1/2,'y');
Sz = sop(1/2,'z');

Mx(1:2,1:numel(p.offsets)) = 0;
My(1:2,1:numel(p.offsets)) = 0;
Mz(1:2,1:numel(p.offsets)) = 0;
for k = 1:2
  
  Amplitude = (Params(k).Flip/Params(k).tp)/(2*pi);
  rho0 = -Sz;
  
  for i = 1:numel(p.offsets)
    
    H = p(k).offsets(i)*Sz + Amplitude*Sy;
    M = -2i*pi*H*Params(k).tp;
    q = sqrt(M(1,1)^2-abs(M(1,2))^2);
    U = cosh(q)*eye(2) + (sinh(q)/q)*M;
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

% Sinc pulse
clearvars -except err
Params.tp = 0.200; % us
Params.Flip = pi;
Params.Type = 'sinc';
Params.zerocross = 0.050;

Opt.Detect = 'Sz';
Opt.Offsets = -100:1:100;
Opt.nBCH = 3;
[t,y,p] = pulse(Params,Opt);

% Calculate inversion profile
Sx = sop(1/2,'x');
Sy = sop(1/2,'y');
Sz = sop(1/2,'z');

Mz(1:numel(p.offsets)) = 0;

rho0 = -Sz;

for i = 1:numel(p.offsets)
  
  rho = rho0;
  
  UPulse = eye(2);
  for j = 1:numel(t)
    H = p.offsets(i)*Sz + real(y(j))*Sx + imag(y(j))*Sy;
    M = -2i*pi*H*(t(2)-t(1));
    q = sqrt(M(1,1)^2 - abs(M(1,2))^2);
    dU = cosh(q)*eye(2) + (sinh(q)/q)*M;
    UPulse = dU*UPulse;
  end
  rho = UPulse*rho*UPulse';
  
  Mz(i) = -2*trace(Sz*rho);
  
end
  
err(2) = ~areequal(Mz,p.Mz,0.5e-1);

% 1st order sech/tanh with frequency offset
clearvars -except err
Params.tp = 0.200; % us
Params.Type = 'sech/tanh';
Params.BW = 100; % MHz
Params.beta = 15;
Params.SweepDirection = -1;
Params.Flip = pi;
Params.TimeStep = 0.0005; % us
Params.CenterFreq = 60; % MHz

t0 = 0:Params.TimeStep:Params.tp;
% Calculate pulse amplitude from flip angle
Qcrit = 5;
BW = Params.BW/tanh(Params.beta/2);
Amplitude = sqrt((Params.beta*BW*Qcrit)/(2*pi*2*Params.tp));
% Amplitude modulation: sech
A = sech((Params.beta/Params.tp)*(t0-Params.tp/2));
% Phase modulation
phi = (Params.BW/(2*tanh(Params.beta/2)))*(Params.tp/Params.beta)*log(cosh((Params.beta/Params.tp)*(t0-Params.tp/2)));
phi = 2*pi*Params.SweepDirection*phi;
% Pulse
y0 = Amplitude*A.*exp(1i*(phi+2*pi*Params.CenterFreq*t0));

[~,~,p] = pulse(Params);

% Inversion profile
Sx = sop(1/2,'x');
Sy = sop(1/2,'y');
Sz = sop(1/2,'z');

rho0 = -Sz;
Mz(1:numel(p.offsets)) = 0;
for i = 1:numel(p.offsets)
  
  H0 = p.offsets(i)*Sz;
  rho = rho0;
  
  UPulse = eye(2);
  for j = 1:numel(t0)
    H = H0 + real(y0(j))*Sx + imag(y0(j))*Sy;   
    M = -2i*pi*H*Params.TimeStep;
    q = sqrt(M(1,1)^2-abs(M(1,2))^2);
    dU = cosh(q)*eye(2) + (sinh(q)/q)*M;
    UPulse = dU*UPulse;
  end 
  rho = UPulse*rho*UPulse';
  
  Mz(i) = -2*trace(Sz*rho);
  
end

err(3) = ~areequal(Mz,p.Mz,1e-3);

err = any(err);

data = [];