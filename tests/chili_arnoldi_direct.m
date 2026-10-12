function ok = test()

% Compare Arnoldi spectrum with direct evaluation for a small damped
% complex-symmetric matrix

rng(1);
N = 60;
H = rand(N); H = H + H.';
G = rand(N); G = G*G.'/N;
A = 1i*H + G;
b = rand(N,1);
b = b/norm(b);
z = -1i*linspace(-2*N,2*N,301);

Opt.Threshold = 1e-10;
[spec,converged] = runprivate('chili_arnoldi',A,b,z,Opt);

I = eye(N);
specref = zeros(size(z));
for k = 1:numel(z)
  specref(k) = b'*((A+z(k)*I)\b);
end

ok = converged && areequal(spec,specref,1e-8,'rel');
