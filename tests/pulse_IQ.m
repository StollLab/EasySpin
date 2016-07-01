function [err,data] = test(opt,olddata)

% Check excitation profiles calculated by pulse()
%--------------------------------------------------------------------------

% Composite pulse 90_0-180_180-270_0
dt = 0.001; % us
tp = 0.008; % us

Params.tp = 6*tp;
v1 = ((pi/2)/tp)/(2*pi); % MHz

I(1:tp/dt) = +1;
I(tp/dt+1:tp/dt+2*tp/dt) = -1;
I(3*tp/dt+1:3*tp/dt+3*tp/dt+1) = +1;

Params.I = v1*I;

% Opt.Offsets = -200:1:200; % MHz
[t,IQ,exprofile] = pulse(Params);

% Calculate inversion profile
Sx = sop(1/2,'x');
Sz = sop(1/2,'z');

Mz(1:numel(exprofile.offsets)) = 0;
rho0 = -Sz;

for i = 1:numel(exprofile.offsets)
  
  H1 = exprofile.offsets(i)*Sz + v1*Sx;
  q1 = sqrt((-2i*pi*tp*H1(1,1))^2-abs((-2i*pi*tp*H1(1,2)))^2);
  U1 = cosh(q1)*eye(2) + (sinh(q1)/q1)*(-2i*pi*tp*H1);
  H2 = exprofile.offsets(i)*Sz - v1*Sx;
  q2 = sqrt((-2i*pi*2*tp*H2(1,1))^2-abs((-2i*pi*2*tp*H2(1,2)))^2);
  U2 = cosh(q2)*eye(2) + (sinh(q2)/q2)*(-2i*pi*2*tp*H2);
  H3 = exprofile.offsets(i)*Sz + v1*Sx;
  q3 = sqrt((-2i*pi*3*tp*H3(1,1))^2-abs((-2i*pi*3*tp*H3(1,2)))^2);
  U3 = cosh(q3)*eye(2) + (sinh(q3)/q3)*(-2i*pi*3*tp*H3);
  U = U3*U2*U1;
  rho = U*rho0*U';
  
  Mz(i) = -2*trace(Sz*rho);
  
end

err(1) = ~areequal(v1*I,real(IQ),1e-12);
err(2) = ~areequal(Mz,exprofile.Mz,1e-12);

err = any(err);

data = [];