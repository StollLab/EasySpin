function ok = test(opt)

% Crash test of perturbation theory with many high-spin nuclei

KerneM.S = 1/2;
KerneM.g = [1.997 1.970 1.965];
KerneM.Nucs = '55Mn, 55Mn, 55Mn, 55Mn';
KerneM.A = [150 160 230; 185 190 265; 220 230 280; 335 326 275];
%noch sind x und z vertauscht - welche Drehung muss man bei der Erzeugung der Orientierungen einfuehren, um das entsprechend auszugleichen?
KerneM.lwEndor = 5;

ExpLsg.Range = [40 250];
ExpLsg.mwFreq = 33.936;
ExpLsg.nPoints = 421;
ExpLsg.ExciteWidth = 200;
OptLsg.Method = 'perturb2';

ExpLsg.Field = 1265;

[rf, spcLsgM] = salt(KerneM, ExpLsg,OptLsg);

if opt.Display
  plot(rf,spcLsgM);
end

ok = true;
