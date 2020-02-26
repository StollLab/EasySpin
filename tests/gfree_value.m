function ok = test(opt)

% Direct value check

a = gfree;
b = 2.00231930436256;
ok = areequal(a,b,1e-14,'abs');
