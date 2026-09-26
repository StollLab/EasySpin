function ok = test()

% Circularly polarized excitation
% (1) crystal: I(circ+) + I(circ-) = 4*I(unpolarized)
% (2) powder: difference spectrum circ+ minus circ- scales with the cosine
%     of the angle between k and B0

Sys.g = [2 2.1 2.2];

Exp.Field = 350;
Exp.SampleFrame = [0 0 0; 0.3 0.7 1.1; 1 2 0.5];
k = [0 sin(0.4) cos(0.4)];
Exp.mwMode = {k 'circular+'};
[~,Ip] = resfreqs_perturb(Sys,Exp);
Exp.mwMode = {k 'circular-'};
[~,Im] = resfreqs_perturb(Sys,Exp);
Exp.mwMode = {k 'unpolarized'};
[~,Iu] = resfreqs_perturb(Sys,Exp);

ok(1) = areequal(Ip+Im,4*Iu,1e-10,'rel');
ok(2) = any(abs(Ip(:)-Im(:))>1e-6*max(Ip(:)));

Sys.lwpp = 20;  % MHz
clear Exp
Exp.Field = 350;
Exp.mwRange = [9.5 11];
Opt.Method = 'perturb';
theta = [0 0.5 1];
for iTheta = numel(theta):-1:1
  k = [0 sin(theta(iTheta)) cos(theta(iTheta))];
  Exp.mwMode = {k 'circular+'};
  [~,spcp] = pepper(Sys,Exp,Opt);
  Exp.mwMode = {k 'circular-'};
  [~,spcm] = pepper(Sys,Exp,Opt);
  dspc(iTheta,:) = spcp-spcm;
end
ratio = dspc(2:end,:)/dspc(1,:);
ok(3) = areequal(ratio(:).',cos(theta(2:end)),1e-6,'rel');
