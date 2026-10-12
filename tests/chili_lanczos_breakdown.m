function ok = test()

% Lanczos detects a serious breakdown and termination, Arnoldi handles the
% breakdown case correctly

Opt.Threshold = 1e-10;
Opt.Lentz = true;
Opt.LanczosBreakdownTol = 1e-5;
z = -1i*linspace(-3,3,101);
directspec = @(A,b) arrayfun(@(z_)b'*((A+z_*eye(size(A)))\b),z);

% Serious breakdown at step 1: y = [0;1;1i] has y.'*y = 0
A = [0.5+2i 1 1i; 1 0.3-1i 0; 1i 0 0.4+0.5i];
b = [1;0;0];
[~,converged,~,~,breakdown] = runprivate('chili_lanczos',A,b,z,Opt);
ok(1) = breakdown && ~converged;
spec = runprivate('chili_arnoldi',A,b,z,Opt);
ok(2) = areequal(spec,directspec(A,b),1e-10,'rel');

% Termination at step 1: b is an eigenvector
A = diag([0.5+2i 0.3-1i 0.4+0.5i]);
[spec,converged,~,nSteps] = runprivate('chili_lanczos',A,b,z,Opt);
ok(3) = converged && nSteps==1 && areequal(spec,directspec(A,b),1e-10,'rel');

% Termination after N steps
rng(1);
N = 5;
H = rand(N); G = rand(N);
A = 1i*(H+H.') + G*G.';
b = rand(N,1); b = b/norm(b);
[spec,converged,~,nSteps] = runprivate('chili_lanczos',A,b,z,Opt);
ok(4) = converged && nSteps==N && areequal(spec,directspec(A,b),1e-8,'rel');
