function ok = test()

% Vinter must have the same shape as the input time axis t, for both
% row and column vectors.

conc = 100;  % µM
lambda = 0.5;

trow = linspace(-1,2,17);
Vrow = dipbackground(trow,conc,lambda);
ok(1) = isequal(size(Vrow),size(trow));

tcol = trow(:);
Vcol = dipbackground(tcol,conc,lambda);
ok(2) = isequal(size(Vcol),size(tcol));

ok(3) = areequal(Vrow(:),Vcol(:),1e-12,'abs');

end
