function ok = test()

value_MHz = 1;
g = 2.005;
value_mT = mhz2mt(value_MHz,g);
value_mT_ref = value_MHz*1e6*planck/bmagn/g/1e-3;

ok = areequal(value_mT,value_mT_ref,1e-14,'rel');

end
