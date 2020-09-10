function ok = test(opt,olddata)

% Check values of Legendre polynomials (M=0)

P{1} = @(z)ones(size(z));
P{2} = @(z)z;
P{3} = @(z)(3*z.^2-1)/2;
P{4} = @(z)(5*z.^3-3*z)/2;
P{5} = @(z)(35*z.^4-30*z.^2+3)/8;
P{6} = @(z)(63*z.^5-70*z.^3+15*z)/8;
P{7} = @(z)(231*z.^6-315*z.^4+105*z.^2-5)/16;
P{8} = @(z)(429*z.^7-693*z.^5+315*z.^3-35*z)/16;
P{9} = @(z)(6435*z.^8-12012*z.^6+6930*z.^4-1260*z.^2+35)/128;

z = 2*rand(1,20)-1;

for k = 1:numel(P)
  L_ = k-1;
  v0 = plegendre(L_,z);
  v1 = P{k}(z);
  ok(k) =  areequal(v0,v1,1e-10,'rel');
end
