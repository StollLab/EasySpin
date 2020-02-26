function ok = test()

% Calling isotopologues() with a spin system without any nuclei

Sys.g = 2;
Iso = isotopologues(Sys);

ok(1) = numel(Iso)==1 && Iso.weight==1;

Iso = isotopologues('',3,[]);

ok(2) = numel(Iso)==1 && Iso.weight==1 && isempty(Iso.n);
