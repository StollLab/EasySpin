function ok = test() %#ok<FNDEF>

% Opt.Binning: delay dimensions computed by binning the peaks in the
% frequency domain, or point by point with fixed delays (exact)
% - binned and point-by-point signals agree within the binning accuracy,
%   for 1D (Carr-Purcell, three incremented intervals), 2D (HYSCORE, both
%   dimensions binned) and mixed choices (one dimension binned)
% - info.binnedDims gives the binned dimensions
% - 'auto': up to 64 points point by point; HYSCORE with many points
%   binned
% - invalid values are rejected: unknown string, wrong number of entries,
%   binning a dimension that is not a delay increment, more than 4
%   incremented intervals, mixed choices with the product rule

Sys.g = gfree;
Sys.Nucs = '1H';
Sys.A = [0.5 1 2];  % MHz
Sys.AFrame = [0 0.6 0];

Exp.mwFreq = 9.8;  % GHz
Exp.Field = unitconvert(Exp.mwFreq,'GHz->mT',gfree);
Exp.MolFrame = [0.3 0.7 0.1];
Exp.CrystalSymmetry = 1;

p90.Flip = pi/2;
p180.Flip = pi;

% Carr-Purcell, delays tau, 2*tau, tau
CP = Exp;
CP.Sequence = {p90 0.1 p180 0.2 p180 0.1};
CP.nPoints = 100;
CP.Dim1 = {'d1' 0.01; 'd2' 0.02; 'd3' 0.01};
CP.Pathways = '-+-';
[~,yStep,info] = saffron(Sys,CP,struct('Binning',false));
ok(1) = isempty(info.binnedDims);
[~,yBin,info] = saffron(Sys,CP,struct('Binning',true));
ok(2) = isequal(info.binnedDims,1);
ok(3) = areequal(yBin,yStep,1e-4,'rel');

% HYSCORE: both dimensions binned, one binned, none binned
HY = Exp;
HY.Sequence = {p90 0.136 p90 0.05 p180 0.05 p90 0.136};
HY.nPoints = [40 40];
HY.Dim1 = {'d2' 0.016};
HY.Dim2 = {'d3' 0.016};
HY.Pathways = '+00-';
[~,yStep,info] = saffron(Sys,HY,struct('Binning',false));
ok(4) = isempty(info.binnedDims);
[~,yBin,info] = saffron(Sys,HY,struct('Binning',true));
ok(5) = isequal(info.binnedDims,[1 2]) && areequal(yBin,yStep,1e-4,'rel');
[~,yMix1,info] = saffron(Sys,HY,struct('Binning',[true false]));
ok(6) = isequal(info.binnedDims,1) && areequal(yMix1,yStep,1e-5,'rel');
[~,yMix2,info] = saffron(Sys,HY,struct('Binning',[false true]));
ok(7) = isequal(info.binnedDims,2) && areequal(yMix2,yStep,1e-5,'rel');

% 'auto': up to 64 points exact, HYSCORE with many points binned
E = CP;
E.nPoints = 64;
[~,~,info] = saffron(Sys,E);
ok(8) = isempty(info.binnedDims);
E = HY;
E.nPoints = [128 128];
[~,~,info] = saffron(Sys,E);
ok(9) = isequal(info.binnedDims,[1 2]);

% invalid values
E = HY;
E.Dim2 = {'p3.Flip' pi/40};  % not a delay increment
Sys2 = Sys;
Sys2.Nucs = '1H,1H';
Sys2.A = [0.5 1 2; 1 1 3];
Sys2.AFrame = [0 0.6 0; 0 0 0];
CP5 = Exp;
CP5.Sequence = {p90 0.1 p180 0.2 p180 0.2 p180 0.2 p180 0.1};
CP5.nPoints = 100;
CP5.Dim1 = {'d1' 0.01; 'd2,d3,d4' 0.02; 'd5' 0.01};
CP5.Pathways = '-+-+-';
cases = {Sys,HY,struct('Binning','fast'),'must be ''auto''';
  Sys,HY,struct('Binning',[true false true]),'one true/false per dimension';
  Sys,E,struct('Binning',[false true]),'cannot be binned';
  Sys,CP5,struct('Binning',true),'limited to 4 incremented intervals';
  Sys2,HY,struct('Binning',[true false],'ProductRule',true),'all dimensions must be binned or all computed point by point'};
for k = 1:size(cases,1)
  try
    saffron(cases{k,1:3});
    ok(end+1) = false; %#ok<AGROW>
  catch err
    ok(end+1) = contains(err.message,cases{k,4}); %#ok<AGROW>
  end
end

end
