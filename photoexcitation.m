tdm_angles = [pi/2 pi/2];
depol = 0;
mode = 'x';
ori = [0 0 pi/2];

photoexcitationweight(tdm_angles,mode,ori,depol)

%%

clear, clc

tdm_angles = [90 90]*pi/180; % Euler angle for tdm vector, in molecular frame
depol_fraction = 0.1; % depolarization contribution, between 0 and 1
mode = 'perp'; % light E-field mode, 'perp' or 'para' relative to B0

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
Triplet.D = [900 160];
Triplet.lwpp = 1;
Triplet.tdm = [0 pi/2];
Triplet.Pop = [1 0 0];

Exp.mwFreq = 9.5;
Exp.Range = [280 400];
Exp.Harmonic = 0;

Opt.GridSymmetry = 'D2h';

Exp.lightPolDir = '';
[B,spc_iso] = pepper(Triplet,Exp,Opt);

Exp.lightPolarization = 0.7;
Exp.lightPolDir = 'x';
[B,spc_perp] = pepper(Triplet,Exp,Opt);
Exp.lightPolDir = 'z';
[B,spc_para] = pepper(Triplet,Exp,Opt);

spc_sum = (spc_para+2*spc_perp)-(1-Exp.lightPolarization)*spc_iso*2;

Exp.lightPolarization = 0;
Exp.lightPolDir = 'x';
[B,spc_unpol] = pepper(Triplet,Exp,Opt);


subplot(2,1,1)
h = plot(B,spc_iso,B,spc_perp,B,spc_para);
h(1).Color = [1 1 1]*0.8;
legend('iso','perp','para')
subplot(2,1,2)
plot(B,spc_iso,B,spc_sum,'--',B,spc_unpol);
legend('iso','para+2*perp','unpol')
