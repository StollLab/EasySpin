function ok = test() %#ok<FNDEF>

% Linear delay increments of any form are simulated by the manifold engine
% (kind-1 dimensions): each interval gets longer by a fixed multiple of the
% time increment per step. The manifold engine, computing each point with
% fixed delays (Opt.Binning = false, exact), agrees with the full
% engine (single orientation, ideal pulses).
% - several incremented delays in one interval: 2p ESEEM with split
%   delays, the same as 2p ESEEM with the summed delays
% - three incremented intervals in one dimension, [1 1 1]
% - two dimensions that increment the same intervals
% - Carr-Purcell with two refocusing pulses, delays tau, 2*tau, tau: steps
%   of different size, [1;2;1]
% - a moving pulse ('p2.Position'): the delay before it gets longer, the
%   one after it shorter
% - a delay that gets shorter (negative step)
% - an interval that changes along both dimensions
% - all echoes are refocused at the detection point at all points, so that
%   the high-field approximation of the manifold engine does not show up
%   in the phase of the signal
% - five incremented intervals (Carr-Purcell with four refocusing pulses):
%   more than binning allows, so Opt.Binning = 'auto' steps it
% - with Opt.Binning = 'auto', Opt.Method = 'auto' uses the manifold
%   engine

Sys.g = gfree;
Sys.Nucs = '1H';
Sys.A = [0.5 1 2];  % MHz
Sys.AFrame = [0 0.6 0];

Exp.mwFreq = 9.8;  % GHz
Exp.Field = unitconvert(Exp.mwFreq,'GHz->mT',gfree);
Exp.MolFrame = [0.3 0.7 0.1];
Exp.CrystalSymmetry = 1;
Exp.DetPhase = -pi/2;

p90.Flip = pi/2;
p180.Flip = pi;

% several incremented delays in one interval
Exps{1} = Exp;
Exps{1}.Sequence = {p90 0.1 0.05 p180 0.1 0.05};
Exps{1}.nPoints = 64;
Exps{1}.Dim1 = {'d1,d2,d3,d4' 0.01};
Exps{1}.Pathways = '+-';

% three incremented intervals in one dimension
Exps{2} = Exp;
Exps{2}.Sequence = {p90 0.12 p90 0.1 p90 0.12};
Exps{2}.nPoints = 80;
Exps{2}.Dim1 = {'d1,d2,d3' 0.01};
Exps{2}.Pathways = '+0-';

% two dimensions incrementing the same intervals
Exps{3} = Exp;
Exps{3}.Sequence = {p90 0.12 p90 0.1 p90 0.12};
Exps{3}.nPoints = [6 40];
Exps{3}.Dim1 = {'d1,d3' 0.04};
Exps{3}.Dim2 = {'d1,d3' 0.01};
Exps{3}.Pathways = '+0-';

% Carr-Purcell, delays tau, 2*tau, tau
Exps{4} = Exp;
Exps{4}.Sequence = {p90 0.1 p180 0.2 p180 0.1};
Exps{4}.nPoints = 60;
Exps{4}.Dim1 = {'d1' 0.01; 'd2' 0.02; 'd3' 0.01};
Exps{4}.Pathways = '-+-';

% moving pulse, with the last delay shortened so that the echo stays
% refocused at the detection point
Exps{5} = Exp;
Exps{5}.Sequence = {p90 0.1 p180 0.3 p180 0.2};
Exps{5}.nPoints = 10;
Exps{5}.Dim1 = {'p2.Position' 0.01; 'd3' -0.02};
Exps{5}.Pathways = '-+-';

% a delay that gets shorter (stimulated echo with shrinking mixing time)
Exps{6} = Exp;
Exps{6}.Sequence = {p90 0.12 p90 0.3 p90 0.12};
Exps{6}.nPoints = 25;
Exps{6}.Dim1 = {'d1,d3' 0.01; 'd2' -0.01};
Exps{6}.Pathways = '+0-';

% an interval that changes along both dimensions
Exps{7} = Exp;
Exps{7}.Sequence = {p90 0.12 p90 0.1 p90 0.12};
Exps{7}.nPoints = [8 10];
Exps{7}.Dim1 = {'d2' 0.02};
Exps{7}.Dim2 = {'d1,d2,d3' 0.01};
Exps{7}.Pathways = '+0-';

% five incremented intervals, more points than the exact limit of 64
Exps{8} = Exp;
Exps{8}.Sequence = {p90 0.1 p180 0.2 p180 0.2 p180 0.2 p180 0.1};
Exps{8}.nPoints = 72;
Exps{8}.Dim1 = {'d1' 0.01; 'd2,d3,d4' 0.02; 'd5' 0.01};
Exps{8}.Pathways = '-+-+-';

tol = 1e-5;  % high-field approximation of the manifold engine
ok = false(numel(Exps),4);
yMans = cell(1,numel(Exps));
for k = 1:numel(Exps)
  E = Exps{k};
  [seq,reasons] = runprivate('sf_parsesequence',E,false);
  ok(k,1) = all(seq.kind==1) && isempty(reasons);
  [~,yMan,info] = saffron(Sys,E,struct('Binning',false));
  yMans{k} = yMan;
  ok(k,2) = strcmp(info.Method,'manifold') && isempty(info.binnedDims);
  E.DetWindow = 0;
  [~,yFull] = saffron(Sys,E,struct('Method','full'));
  ok(k,3) = areequal(yMan,yFull,tol,'rel');
  % the signal is modulated much more than the tolerance
  ok(k,4) = std(real(yMan(:)))/max(abs(yMan(:))) > 10*tol;
end
ok = ok(:);

% Split delays vs 2p ESEEM with the summed delays
E = Exp;
E.Sequence = {p90 0.15 p180 0.15};
E.nPoints = 64;
E.Dim1 = {'d1,d2' 0.02};
E.Pathways = '+-';
[~,y2p] = saffron(Sys,E);
ok(end+1) = areequal(yMans{1},y2p,1e-12,'rel');

% Opt.Method = 'auto' and Opt.Binning = 'auto': manifold engine; five
% incremented intervals are stepped, not binned
[~,y5,info] = saffron(Sys,Exps{8});
ok(end+1) = strcmp(info.Method,'manifold') && isempty(info.binnedDims);
ok(end+1) = areequal(y5,yMans{8},1e-12,'rel');

end
