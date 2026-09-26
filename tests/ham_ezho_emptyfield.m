function ok = test()

% Assert that ham_ezho with an empty field returns the tensors.

rng(5);

Sys.S = 1;
Sys.Ham112 = rand(1,5);

cG = ham_ezho(Sys,[]);

ok = iscell(cG) && numel(cG)==2 && numel(cG{2})==3;
