function ok = test()

amu_val = amu;
amu_ref = 1.66053906892e-27;

ok = areequal(amu_val,amu_ref,1e-12,'rel');
