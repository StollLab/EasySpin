function ok = test()

% Test correct sign of emission/absorption for non-Boltzmann populations
% Low-field line must be in emission, the high-field line in absorption

Sys.S = 1;
Sys.g = [2 2 2];
D = 0.06; E = 0.001;
Sys.D = 100*clight/1e6*[-D/3 + E, -D/3 - E, 2*D/3];

Exp = struct('mwFreq',9.67739,'nPoints',426,'Range',[260 430]);
Exp.Harmonic = 0;
pop = [1 0 0];
Sys.initState = {pop,'zerofield'};

Exp.CrystalOrientation = [0 0 0];
[Bpara,Apara] = resfields(Sys,Exp);
[Bpara,idx] = sort(Bpara);
Apara = Apara(idx);

Exp.CrystalOrientation = [0 pi/2 0];
[Bperp,Aperp] = resfields(Sys,Exp);
[Bperp,idx] = sort(Bperp);
Aperp = Aperp(idx);

if numel(Bperp)==2 && numel(Bpara)==2
  ok = Apara(1)>0 && Apara(2)<0 && Aperp(1)>0 && Aperp(2)<0;
else
  ok = false;
end
