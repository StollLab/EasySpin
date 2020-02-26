function ok = test()

% Syntax test
q = 1000;
z = rand(1,q);
y = plegendre(16,-5,z);
plegendre(16,-5,z);
plegendre(3,z);

ok = true;
