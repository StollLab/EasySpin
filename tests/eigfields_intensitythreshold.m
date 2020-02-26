function ok = test()

%===============================================================================
% eigfields should use Opt.Threshold to screen transitions
% no matter whether the user requests intensities or not.
%===============================================================================

Exp.mwFreq=9.65;

Sys.g = gfree;
Sys.Nucs = '31P,31P';
Sys.A = mt2mhz([10 20]/10);

Op.Threshold = 0.001;

[Pos1,Int] = eigfields(Sys,Exp,Op);
Pos2 = eigfields(Sys,Exp,Op);

ok = areequal(Pos1,Pos2,1e-10,'abs');

