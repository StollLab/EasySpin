function ok = test()

value_mT = 1;
g = 2.005;
value_MHz = mt2mhz(value_mT,g);
value_MHz_ref = value_mT*1e-3*bmagn*g/1e6/planck;

ok = areequal(value_MHz,value_MHz_ref,1e-14,'rel');

end
