function [err,data] = test(opt,olddata)

% Assure that running esfit with multiple, sequential algorithms is 
% successful and yields a good fit for a custom function that generates a 
% 3D array of data.

FitOpt.nTrials = 200;

% number of grid points in each direction
N = 50;

rng_(1)

err = false;

FitOpt.PrintLevel = 0;

AddNoise = 1;

% generate Euler angle grids
agrid = linspace(-180, 180, N)/180*pi;
bgrid = linspace(-90, 90, N)/180*pi;
ggrid = linspace(-180, 180, N)/180*pi;
[Agrid, Bgrid, Ggrid] = meshgrid(agrid, bgrid, ggrid);
grids = {Agrid, Bgrid, Ggrid};

% number of Gaussians
nComps = 1;

% generate parameter and fit range values
for iComp = 1:nComps
  x(iComp,1:2) = [45, 45]/180*pi;  % alpha: alpha0, fwhm
  x(iComp,3:4) = [45, 45]/180*pi;  % beta: beta0, fwhm
  x(iComp,5:6) = [45, 45]/180*pi;  % beta: beta0, fwhm
  
  % alpha-direction
  Vary{iComp,1} = [0, 90]/180*pi;  % alpha0
  Vary{iComp,2} = [10, 60]/180*pi;  % fwhm
  
  % beta-direction
  Vary{iComp,3} = [0, 90]/180*pi;  % beta0
  Vary{iComp,4} = [10, 60]/180*pi;  % fwhm
  
  % gamma-direction
  Vary{iComp,5} = [0, 90]/180*pi;  % beta0
  Vary{iComp,6} = [10, 60]/180*pi;  % fwhm
end

% data to fit
data = gaussfunc(grids, x);

if AddNoise
  noiseamp = 0.1*max(data(:));
  noise = noiseamp*randn(N,N,N);
  data = data + noise;
end

% initial guess
x0 = [55, 55,...
      55, 55,...
      55, 55]/180*pi;

% transform to [-1,1] to use esfit algorithm
x0 = invtransformVars(x0, Vary);

nParameters = numel(x0);

xfit = esfit_montecarlo(@getresid, nParameters, FitOpt, data, grids, Vary);  % info is not set to a value in esfit_montecarlo

[xfit, ignore] = esfit_levmar(@getresid, xfit, [], data, grids, Vary);

% transform back for plotting
xfit = transformVars(xfit, Vary);

sd = (data - gaussfunc(grids, xfit)).^2;
rmsd = real(sqrt(mean(sd(:))));

err = rmsd>0.5;

data = [];

end

% helper functions
% -------------------------------------------------------------------------

function out = gaussfunc(grids, params)

X = grids{1};
Y = grids{2};
Z = grids{3};

out = 0;

for iComp=1:size(params,1)
  X0 = params(iComp, 1);
  Xfwhm = params(iComp, 2);
  
  Y0 = params(iComp, 3);
  Yfwhm = params(iComp, 4);
  
  Z0 = params(iComp, 5);
  Zfwhm = params(iComp, 6);
  
  out = out + gaussian(X, X0, Xfwhm).*...
              gaussian(Y, Y0, Yfwhm).*...
              gaussian(Z, Z0, Zfwhm);
end

% int = sum(sum(out,1),2)*(X(1,2)-X(1,1))*(Y(2,1)-Y(1,1));
% out = out/int;

end

function [rmsd, resid, simdata] = getresid(x, zdata, grids, Vary)

global smallestError

xp = transformVars(x, Vary);

simdata = gaussfunc(grids, xp);
resid = zdata - simdata;
rmsd = real(sqrt(mean(mean(mean(resid.^2,1),2),3)));

isNewBest = rmsd < smallestError;

if isNewBest
  smallestError = rmsd;
end

end

function xp = invtransformVars(x, Vary)
% transform variables from Vary domains to [-1,1] domain

nComps = size(Vary, 1);

VaryCents = zeros(nComps, 6);  % center of each Vary domain
VaryRanges = zeros(nComps, 6);  % range of each Vary domain

for iComp=1:nComps
  VaryCents(iComp, :) = [mean(Vary{iComp,1}), mean(Vary{iComp,2}), ...
                         mean(Vary{iComp,3}), mean(Vary{iComp,4}), ...
                         mean(Vary{iComp,5}), mean(Vary{iComp,6})];
  VaryRanges(iComp, :) = [findRange(Vary{iComp,1}), ...
                          findRange(Vary{iComp,2}), ...
                          findRange(Vary{iComp,3}), ...
                          findRange(Vary{iComp,4}), ...
                          findRange(Vary{iComp,5}), ...
                          findRange(Vary{iComp,6})];
end

xp = (x - VaryCents)./VaryRanges;

end

function xp = transformVars(x, Vary)
% transform variables from [-1,1] domain to Vary domains

nComps = size(Vary, 1);
x = reshape(x, [nComps,6]);  % esfit algorithms set x = x(:), so undo this 
                             % here

VaryCents = zeros(nComps, 6);  % center of each Vary domain
VaryRanges = zeros(nComps, 6);  % range of each Vary domain

for iComp=1:nComps
  VaryCents(iComp, :) = [mean(Vary{iComp,1}), mean(Vary{iComp,2}), ...
                         mean(Vary{iComp,3}), mean(Vary{iComp,4}), ...
                         mean(Vary{iComp,5}), mean(Vary{iComp,6})];
  VaryRanges(iComp, :) = [findRange(Vary{iComp,1}), ...
                          findRange(Vary{iComp,2}), ...
                          findRange(Vary{iComp,3}), ...
                          findRange(Vary{iComp,4}), ...
                          findRange(Vary{iComp,5}), ...
                          findRange(Vary{iComp,6})];
end

xp = x.*VaryRanges + VaryCents;

end

function range = findRange(x)

range = (max(x) - min(x))/2;

end