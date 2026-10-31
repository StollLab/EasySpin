function ok = test()

% Matrix T1 and T2 for a spin-1/2: T1(i,j) gives the rate 1/(2*T1(i,j))
% from state j to state i, lower triangle of T2 is ignored

H = diag([1 -1]);  % MHz
tol = 1e-12;

% Different rates up and down
Sys.Spins = 1/2;
Sys.nElectrons = 1;
Sys.T1 = [0 1.5; 1 0];
Sys.T2 = 0.5;

w12 = 1/(2*1.5);  % from state 2 to state 1
w21 = 1/(2*1);    % from state 1 to state 2
Ref = zeros(4);
Ref([1 4],[1 4]) = [w21 -w12; -w21 w12];
Ref(2,2) = 1/Sys.T2;
Ref(3,3) = 1/Sys.T2;

Gamma = runprivate('s_relaxationsuperoperator',Sys,H);
ok(1) = areequal(Gamma,Ref,tol,'abs');

% Only upper triangle: same rate in both directions
Sys.T1 = [0 1.5; 0 0];
Sys.T2 = [0 0.5; 0.3 0];  % 0.3 is ignored

w = 1/(2*1.5);
Ref = zeros(4);
Ref([1 4],[1 4]) = [w -w; -w w];
Ref(2,2) = 1/0.5;
Ref(3,3) = 1/0.5;

Gamma = runprivate('s_relaxationsuperoperator',Sys,H);
ok(2) = areequal(Gamma,Ref,tol,'abs');

% No T1
Sys.T1 = zeros(2);
Ref([1 4],[1 4]) = 0;
Gamma = runprivate('s_relaxationsuperoperator',Sys,H);
ok(3) = areequal(Gamma,Ref,tol,'abs');
