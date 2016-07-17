function [err,data] = test(opt,olddata)

% Test all g input expansions and transformations
%=================================================

% isotropic, one electron
Sys.S = 1/2;
Sys.g = 2;

Q = runprivate('validatespinsys',Sys);
err(1) = ~areequal(Q.g,[2 2 2]);

% isotropic, several electrons
Sys.S = [1 1];
Sys.g = [2 3];
Sys.ee = [0 0 0];
Q = runprivate('validatespinsys',Sys);
err(2) = ~areequal(Q.g,[2 2 2; 3 3 3]);

Sys.S = [1 1];
Sys.g = [2;3];
Sys.ee = [0 0 0];
Q = runprivate('validatespinsys',Sys);
err(3) = ~areequal(Q.g,[2 2 2; 3 3 3]);

% anisotropic, one electron
clear Sys
Sys.g = [2 3 4];
Q = runprivate('validatespinsys',Sys);
err(4) = ~areequal(Q.g,Sys.g);

% anisotropic, several electron
clear Sys
Sys.S = [1/2,1/2];
Sys.g = [2 3 4; 5 6 7];
Sys.ee = [1 1 1];
Q = runprivate('validatespinsys',Sys);
err(5) = ~areequal(Q.g,Sys.g);

% absent, one electron
clear Sys
Sys.S = 1;
Q = runprivate('validatespinsys',Sys);
err(6) = ~areequal(Q.g,[1 1 1]*gfree);

% absent, several electron
clear Sys
Sys.S = [1/2,3/2];
Sys.ee = [0 0 0];
Q = runprivate('validatespinsys',Sys);
err(7) = ~areequal(Q.g,[1 1 1; 1 1 1]*gfree);

data = [];
