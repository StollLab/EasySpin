function ok = test()

g1H = gamman('1H');
ref = nucgval('1H')*nmagn/hbar;

ok = areequal(g1H,ref,1e-10,'rel');

end
