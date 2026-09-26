function ok = test()

% Assert that ham with higher-order Zeeman terms works with 'sparse', without
% a field, and with only zero-field terms

rng(5);

Sys1.S = 1;
Sys1.Ham112 = rand(1,5);
[H0,mux,muy,muz] = ham(Sys1,[],'sparse');
ok(1) = issparse(H0) && issparse(mux) && issparse(muy) && issparse(muz);

Sys0.S = 1;
Sys0.Ham022 = rand(1,5);
[H0,mux,muy,muz] = ham(Sys0);
ok(2) = areequal(H0,ham_ezho(Sys0,[0 0 0]),1e-10,'abs') && ~any(mux(:));
