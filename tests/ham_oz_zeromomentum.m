function ok = test()

% Make sure all-zero matrix is returned if the orbital angular momemtum is 0.

Sys.S = 1/2;
Sys.L = 0;
Sys.soc = 1;
H = ham_oz(Sys,[0 0 1000]);

ok = all(H(:)==0) && size(H,1)==2;
