function ok = test()

%===============================================================================
% resfields_eig should use Opt.Threshold to screen transitions
% no matter whether the user requests intensities or not.
%===============================================================================

Exp.mwFreq=9.65;

Sys.g = gfree;
Sys.Nucs = '31P,31P';
Sys.A = unitconvert([10 20]/10,'mT->MHz');

Op.Threshold = 0.001;

[Pos1,Int] = resfields_eig(Sys,Exp,Op);
Pos2 = resfields_eig(Sys,Exp,Op);

ok = areequal(Pos1,Pos2,1e-10,'abs');

