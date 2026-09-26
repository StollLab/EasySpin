function ok = test()

% Assert that ham_ezho returns sparse matrices only if requested, also if
% there are no terms of the requested order in B.

rng(5);

Sys.S = 1;
Sys.Ham112 = rand(1,5);
B = [1 2 3];

ok(1) = ~issparse(ham_ezho(Sys,B));
ok(2) = issparse(ham_ezho(Sys,B,[],'sparse'));
ok(3) = ~issparse(ham_ezho(Sys,B,[],'',0));  % no terms of order 0
ok(4) = issparse(ham_ezho(Sys,B,[],'sparse',0));
