clear, clc

tdm_angles = [90 90]*pi/180; % Euler angle for tdm vector, in molecular frame
depol_fraction = 0.1; % depolarization contribution, between 0 and 1
mode = 'perp'; % laser E-field mode, 'perp' or 'para' relative to B0

% Euler angles of lab frame, relative to mol frame
phi = 0;
theta = 0;
chi = 0;
chilist = linspace(0,pi,101)
for ichi = 1:numel(chilist)
  w(ichi) = photoexcitationweight(tdm_angles,mode,[phi theta chilist(ichi)],depol_fraction)
end

%%

clear, clc

Triplet.S = 1;
Triplet.D = 400;
Triplet.lwpp = 0.5;
Triplet.tdm = [0 0];

Exp.mwFreq = 9.5;
Exp.Range = [320 360];
Exp.Harmonic = 0;

Exp.laserPolDir = '';
[B,spc] = pepper(Triplet,Exp);

Exp.laserDepolarization = 0;
Exp.laserPolDir = 'x';
[B,spc_perp] = pepper(Triplet,Exp);
Exp.laserPolDir = 'z';
[B,spc_para] = pepper(Triplet,Exp);

spc_iso = (spc_para+2*spc_perp)*(1-Exp.laserDepolarization);

subplot(2,1,1)
plot(B,spc,B,spc_perp,B,spc_para);
legend('-','perp','para')
subplot(2,1,2)
plot(B,spc,B,spc_iso,'--');
legend('-','sum')
