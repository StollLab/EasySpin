function [err,data] = test(opt,olddata)

% Calling isotopologues() with a spin system without any nuclei
%-------------------------------------------------------------------------------

Sys.g = 2;
Iso = isotopologues(Sys);

err = numel(Iso)~=1 || Iso.weight~=1;

Iso = isotopologues('',3,[]);

err = err || numel(Iso)~=1 || Iso.weight~=1 || ~isempty(Iso.n);

data = [];
