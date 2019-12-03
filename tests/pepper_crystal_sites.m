function [err,data] = test(opt,olddata)

% Simulate individual crystal sites

Sys.g = [2 2.1 2.2];
Sys.lwpp = 1;

Exp.mwFreq = 9.5;
Exp.Range = [200 400];

Exp.CrystalSymmetry = 34;
Exp.CrystalOrientation = rand(1,3)*pi;

Opt.Sites = []; % all sites
[x,y0] = pepper(Sys,Exp);

y1 = 0;
nSites = 4;
for s=1:nSites
  Opt.Sites = 1;
  y1 = y1 + pepper(Sys,Exp);
end
y1 = y1/nSites;

if opt.Display
  plot(x,y0,x,y1);
  legend('all sites at once','each site separately');
end

err = ~areequal(y0,y1,1e-10,'rel');
data = [];
