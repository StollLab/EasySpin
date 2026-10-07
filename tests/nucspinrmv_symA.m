function ok = test()

% Removing nuclei with symmetric A matrices [xx yy zz xy xz yz]

A = [1 2 3 0.1 0.2 0.3; 4 5 6 0.4 0.5 0.6; 7 8 9 0.7 0.8 0.9];
Sys = struct('S',1/2,'Nucs','1H,15N,13C','A',A);

Sys1 = nucspinrmv(Sys,2);
ok(1) = isequal(Sys1.A,A([1 3],:)) && strcmp(Sys1.Nucs,'1H,13C');

Sys2 = nucspinrmv(Sys,[1 3]);
ok(2) = isequal(Sys2.A,A(2,:)) && strcmp(Sys2.Nucs,'15N');
