function ok = test(opt)

% A powder spectrum of an spin system with axial g tensor should be
% independent of the region of orientational integration.

Sys.g = [1.9,2.3]; % axial g tensor
Sys.lw = 1;

Exp.mwFreq = 9.5;
Exp.Range = [285 365];
Exp.Harmonic = 0;

Opt.Method = 'perturb1';
Opt.GridSize = [19 5];

Symmetry = {'Dinfh','D6h','D4h','Oh','D3d','Th','D2h',...
    'C4h','C6h','S6','C2h','Ci','C1'};

for k = 1:numel(Symmetry)
  Opt.GridSymmetry = Symmetry{k};
  [x,y(k,:)] = pepper(Sys,Exp,Opt);
end

yref = y(1,:);   % spectrum with Dinfh grid is reference

scale = max(yref);
yref = yref/scale;
y = y/scale;

for k = 1:numel(Symmetry)
  ok(k) = max(abs(y(k,:)-yref))<0.01;
end

if opt.Display
  plot(x,y);
  title('Axial symmetry invariance');
  xlabel('magnetic field [mT]');
  ylabel('normalized amplitude');
  legend(Symmetry{:});
end

