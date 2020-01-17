function err = test(opt)

% Direct value check

a = gfree;
b = 2.00231930436256;
err = areequal(a,b,1e-14,'abs');
