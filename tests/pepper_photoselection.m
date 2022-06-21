function ok = test(opt)

% Test photoselection. A spectrum with isotropic light excitation
% (experimentally not possible) should be equal to the weighted sum of 
% spectra with excitation with E-field perpendicular and parallel to B0.

Triplet.S = 1;
Triplet.D = [900 160];  % MHz
Triplet.lwpp = 1;  % mT
Triplet.Pop = [1 0 1];
Triplet.tdm = 'y';  % transition dipole moment along yMol

Exp.mwFreq = 9.5;  % GHz
Exp.Range = [280 400];  % mT
Exp.Harmonic = 0;

Opt.GridSymmetry = 'D2h';

iso = 0.2;
Exp.lightScatter = iso;  % isotropic contribution

Exp.lightMode = '';
[B,spc_iso] = pepper(Triplet,Exp,Opt);
Exp.lightMode = 'perp';
[B,spc_perp] = pepper(Triplet,Exp,Opt);
Exp.lightMode = 'para';
[B,spc_para] = pepper(Triplet,Exp,Opt);

spc_sum = (spc_para+2*spc_perp) - 2*iso*spc_iso;

if opt.Display
  plot(B,spc_iso,B,spc_sum)
end

ok = areequal(spc_iso,spc_sum,1e-10,'rel');

end
