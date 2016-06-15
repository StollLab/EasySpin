function [err,data] = test(opt,olddata)

%====================================================
% orisel with C2h symmetry
%====================================================
% should compute as many knots as salt will use.

% Definition of the spin system
System = struct('S',1/2,'g',[2.1 2.111 2.219],'Nucs','13C','A',[0.3 19.68 24.2],'AFrame',[0 -20 0]);

% Definition of the experiment
Experiment = struct('Range',[-1 1]*5+13.3,'mwFreq',9.5,'ExciteWidth',50);
Experiment.Field = 310;

nKnots = 5;
Symmetry = 'C2h';

% ENDOR simulation options
Options = struct('Verbosity',0,'nKnots',nKnots,'Symmetry',Symmetry);
Options.OriWeights = orisel(System,Experiment,Options);

[phi,theta] = sphgrid(Symmetry,nKnots);

data = [];
err = numel(Options.OriWeights)~=numel(phi);

if opt.Display
  fprintf('  Number of orientations from orisel: %d\n  Number from sphgrid: %d\n',...
    numel(Options.OriWeights), numel(phi));
end

if isempty(err)
  [x,y] = salt(System,Experiment,Options);
end
