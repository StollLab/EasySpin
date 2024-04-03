function ok = test()

% Check whether resfreqs_matrix returns the microwave frequency if
% supplied with resonance fields determined by resfreqs.

Sys.S = 1;
Sys.D = 30e3*0.2*[1 0.2];
Sys.Nucs = '1H';
Sys.A = [100 150 -200];

Exp1.mwFreq = 10;
Exp1.Range = [150 400];

Exp1.SampleFrame = [0.7 0.43 0.1231]*pi;
Exp2.SampleFrame = Exp1.SampleFrame;
[B,~,~,Tr] = resfields(Sys,Exp1);
for iB = 1:numel(B)
  Exp2.Field = B(iB);
  Opt.Transitions = Tr(iB,:);
  nu(iB) = resfreqs_matrix(Sys,Exp2,Opt)/1e3;
end
ok = all(nu-Exp1.mwFreq<1e-6);
