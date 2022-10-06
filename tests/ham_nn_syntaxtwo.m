function ok = test()

% Syntax test, 2 nuclei

Sys.Nucs = '1H,13C';
Sys.A = [10 3];

% Isotropic coupling
Sys.nn = 0.1;

ham_nn(Sys);
ham_nn(Sys,[1 2]);
H = ham_nn(Sys);
H = ham_nn(Sys,[1 2]);

% Principal values
Sys.nn = 0.1*[1 1 1];

ham_nn(Sys);
ham_nn(Sys,[1 2]);
H = ham_nn(Sys);
H = ham_nn(Sys,[1 2]);

% Full tensor
Sys.nn = 0.1*eye(3);

ham_nn(Sys);
ham_nn(Sys,[1 2]);
H = ham_nn(Sys);
H = ham_nn(Sys,[1 2]);

ok = true;
