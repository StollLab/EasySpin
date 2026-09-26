function ok = test()

% Assert that a term with lS=0 throws an error.

rng(5);

Sys.S = 1;
Sys.Ham202 = rand(1,5);

try
  ham_ezho(Sys,[0 0 1]);
  ok = false;
catch err
  ok = contains(err.message,'lS must be');
end
