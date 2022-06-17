clc
rng(4)
tdmOri = rand(1,2)*pi;
ori = [0 0 pi/2];
kOri = 'y';
alpha = 0;

photoselect(tdmOri,ori,kOri,alpha)

%%

clear, clc

tdmOri = 'y';
k = [pi/2 pi/2];
alpha = pi/2;

% Euler angles of lab frame, relative to mol frame
phi = 0;
theta = 0;
chilist = linspace(0,pi,101);
for ichi = 1:numel(chilist)
  ori = [phi theta chilist(ichi)];
  w(ichi) = photoselect(tdmOri,ori,k,alpha);
end

plot(chilist,w)

%%

clear, clc

Triplet.S = 1;
Triplet.D = [900 160];
Triplet.lwpp = 1;
Triplet.tdm = 'y';
Triplet.Pop = [1 0 0];

Exp.mwFreq = 9.5;
Exp.Range = [280 400];
Exp.Harmonic = 0;

Opt.GridSymmetry = 'D2h';

Exp.lightPolAngle = NaN;
[B,spc_iso] = pepper(Triplet,Exp,Opt);

Exp.lightScatter = 0.3;
Exp.lightPolAngle = pi/2;
[B,spc_perp] = pepper(Triplet,Exp,Opt);
Exp.lightPolAngle = 0;
[B,spc_para] = pepper(Triplet,Exp,Opt);

spc_sum = (spc_para+2*spc_perp)-Exp.lightScatter*spc_iso*2;

Exp.lightPolAngle = NaN;
[B,spc_unpol] = pepper(Triplet,Exp,Opt);
spc_unpol = rescaledata(spc_unpol,spc_sum,'lsq');

subplot(2,1,1)
h = plot(B,spc_iso,B,spc_perp,B,spc_para);
h(1).Color = [1 1 1]*0.8;
legend('iso','perp','para')
subplot(2,1,2)
plot(B,spc_iso,B,spc_sum,'--',B,spc_unpol);
legend('iso','para+2*perp','unpol')
