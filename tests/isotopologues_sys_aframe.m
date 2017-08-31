function [err,data] = test(opt,olddata)

% Isotopologues with AFrame
%-------------------------------------------------------------------------------
Sys.Nucs = 'C';
Sys.n = 3;
Sys.A = [10 10 20];
Sys.AFrame = [1 2 2]*pi/3;

Iso = isotopologues(Sys,0);

err = numel(Iso(1).AFrame)~=0 || ...
      numel(Iso(2).AFrame)~=3 || ...
      numel(Iso(3).AFrame)~=3 || ...
      numel(Iso(4).AFrame)~=3;

err = err || numel(Iso)~=4;

data = [];
