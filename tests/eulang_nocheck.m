function ok = test()

% With nocheck = true, eulang skips all input checks.

rng(5231);

% Valid rotation matrix: same result with and without checks
R = erot(rand(1,3)*pi.*[2 1 2]);
ok(1) = areequal(eulang(R,true),eulang(R),1e-14,'abs');

% Slightly non-orthogonal matrix: with checks, eulang prints a message and
% orthogonalizes; with nocheck, it prints nothing and uses the matrix as is
R = erot([10 40 80]*pi/180) + 0.001*rand(3,3); %#ok<NASGU> used in evalc
out_check = evalc('ang_check = eulang(R);');
out_nocheck = evalc('ang_nocheck = eulang(R,true);');
ok(2) = ~isempty(out_check);
ok(3) = isempty(out_nocheck);
ok(4) = ~areequal(ang_nocheck,ang_check,1e-6,'abs');
