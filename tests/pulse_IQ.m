function [err,data] = test(opt,olddata)

% Check excitation profiles calculated by pulse()
%--------------------------------------------------------------------------

% Composite pulse 90_0-180_180-270_0
dt = 0.001; % us
tp = 0.008; % us

Exp.tp = 6*tp;
v1 = ((pi/2)/tp)/(2*pi); % MHz

I(1:tp/dt) = +1;
I(tp/dt+1:tp/dt+2*tp/dt) = -1;
I(3*tp/dt+1:3*tp/dt+3*tp/dt+1) = +1;

Exp.PulseShape.I = v1*I;

Opt.Offsets = -200:1:200;
[~,y,p] = pulse(Exp,Opt);

% Calculate inversion profile
Sx = sop(1/2,'x');
Sz = sop(1/2,'z');

Mz(1:numel(p.offsets)) = 0;
rho0 = -Sz;

for i = 1:numel(p.offsets)
  
  H1 = p.offsets(i)*Sz + v1*Sx;
  U1 = expm(-2i*pi*H1*tp);
  H2 = p.offsets(i)*Sz - v1*Sx;
  U2 = expm(-2i*pi*H2*2*tp);
  H3 = p.offsets(i)*Sz + v1*Sx;
  U3 = expm(-2i*pi*H3*3*tp);
  U = U3*U2*U1;
  rho = U*rho0*U';
  
  Mz(i) = -2*trace(Sz*rho);
  
end

err(1) = ~areequal(v1*I,real(y),1e-12);
err(2) = ~areequal(Mz,p.Mz,1e-12);

err = any(err);

data = [];