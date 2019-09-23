function [err,data] = test(opt,olddata)

% Check excitation profiles calculated by exciteprofile()
%--------------------------------------------------------------------------

% Rectangular pulses
Params(1).Flip = pi/2;
Params(1).tp = 0.016; % us
Params(1).Phase = pi/2;
Params(2).Flip = pi;
Params(2).tp = 0.032; % us
Params(2).Phase = pi/2;

[t1,IQ1] = pulse(Params(1));
[t2,IQ2] = pulse(Params(2));

offsets = -70:0.5:70; % MHz
[offsets,M1] = exciteprofile(t1,IQ1,offsets);
[offsets,M2] = exciteprofile(t2,IQ2,offsets);

% Calculate inversion profile
Sx = sop(1/2,'x');
Sy = sop(1/2,'y');
Sz = sop(1/2,'z');

Mx(1:2,1:numel(offsets)) = 0;
My(1:2,1:numel(offsets)) = 0;
Mz(1:2,1:numel(offsets)) = 0;
for k = 1:2
  
  Amplitude = (Params(k).Flip/Params(k).tp)/(2*pi);
  rho0 = -Sz;
  
  for i = 1:numel(offsets)
    
    H = offsets(i)*Sz + Amplitude*Sy;
    M = -2i*pi*H*Params(k).tp;
    q = sqrt(M(1,1)^2-abs(M(1,2))^2);
    U = cosh(q)*eye(2) + (sinh(q)/q)*M;
    rho = U*rho0*U';
    
    Mx(k,i) = -2*trace(Sx*rho);
    My(k,i) = -2*trace(Sy*rho);
    Mz(k,i) = -2*trace(Sz*rho);
    
  end
  
end

thr = 1e-12;
suberr(1) = ~areequal(Mx(1,:),M1(1,:),thr,'rel');
suberr(2) = ~areequal(Mx(2,:),M2(1,:),thr,'rel');
suberr(3) = ~areequal(My(1,:),M1(2,:),thr,'rel');
suberr(4) = ~areequal(My(2,:),M2(2,:),thr,'rel');
suberr(5) = ~areequal(Mz(1,:),M1(3,:),thr,'rel');
suberr(6) = ~areequal(Mz(2,:),M2(3,:),thr,'rel');

err(1) = any(suberr);

% Sinc pulse
clear Params Mx My Mz
Params.tp = 0.200; % us
Params.Flip = pi;
Params.Type = 'sinc';
Params.zerocross = 0.050;

Opt.Detect = 'Sz';
Opt.nBCH = 3;
[t,IQ] = pulse(Params,Opt);

offsets = -100:1:100;
[offsets,Mag] = exciteprofile(t,IQ,offsets);

% Calculate inversion profile
Sx = sop(1/2,'x');
Sy = sop(1/2,'y');
Sz = sop(1/2,'z');

Mz(1:numel(offsets)) = 0;

rho0 = -Sz;

for i = 1:numel(offsets)
  
  rho = rho0;
  
  UPulse = eye(2);
  for j = 1:numel(t)
    H = offsets(i)*Sz + real(IQ(j))*Sx + imag(IQ(j))*Sy;
    M = -2i*pi*H*(t(2)-t(1));
    q = sqrt(M(1,1)^2 - abs(M(1,2))^2);
    dU = cosh(q)*eye(2) + (sinh(q)/q)*M;
    UPulse = dU*UPulse;
  end
  rho = UPulse*rho*UPulse';
  
  Mz(i) = -2*trace(Sz*rho);
  
end
  
err(2) = ~areequal(Mz,Mag(3,:),5e-2,'rel');

% 1st order sech/tanh with frequency offset
clear Params Mx My Mz
Params.tp = 0.200; % us
Params.Type = 'sech/tanh';
Params.Frequency = [50 -50] + 60; % MHz
Params.beta = 15;
Params.Flip = pi;
Params.TimeStep = 0.0005; % us

t0 = 0:Params.TimeStep:Params.tp;
% Calculate pulse amplitude from flip angle
Qcrit = 5;
BW = (Params.Frequency(2)-Params.Frequency(1))/tanh(Params.beta/2);
Amplitude = sqrt((Params.beta*BW*Qcrit)/(2*pi*2*Params.tp));
% Amplitude modulation: sech
A = sech((Params.beta/Params.tp)*(t0-Params.tp/2));
% Phase modulation
phi = ((Params.Frequency(2)-Params.Frequency(1))/(2*tanh(Params.beta/2)))*(Params.tp/Params.beta)*log(cosh((Params.beta/Params.tp)*(t0-Params.tp/2)));
phi = 2*pi*phi;
% Pulse
IQ0 = Amplitude*A.*exp(1i*(phi+2*pi*mean(Params.Frequency)*t0));

[t,IQ] = pulse(Params);
[offsets,Mag] = exciteprofile(t,IQ);

% Inversion profile
Sx = sop(1/2,'x');
Sy = sop(1/2,'y');
Sz = sop(1/2,'z');

rho0 = -Sz;
Mz(1:numel(offsets)) = 0;
for i = 1:numel(offsets)
  
  H0 = offsets(i)*Sz;
  rho = rho0;
  
  UPulse = eye(2);
  for j = 1:numel(t0)
    H = H0 + real(IQ0(j))*Sx + imag(IQ0(j))*Sy;   
    M = -2i*pi*H*Params.TimeStep;
    q = sqrt(M(1,1)^2-abs(M(1,2))^2);
    dU = cosh(q)*eye(2) + (sinh(q)/q)*M;
    UPulse = dU*UPulse;
  end 
  rho = UPulse*rho*UPulse';
  
  Mz(i) = -2*trace(Sz*rho);
  
end

err(3) = ~areequal(Mz,Mag(3,:),1e-3,'rel');

err = any(err);

data = [];