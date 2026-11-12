function ok = test() %#ok<FNDEF>

% Linear delay increments that do not fit a binning scheme of sf_peaks are
% computed point by point in the subspace engine, instead of switching to
% the full-space engine
% - several incremented delays in one interval: 2p ESEEM with split delays,
%   the same as 2p ESEEM with the summed delays
% - three incremented intervals in one dimension (scheme [1 1 1])
% - two dimensions that increment the same interval: the larger one stays
%   a delay increment (kind 1), the other is computed point by point
% - Opt.Method = 'auto' uses the subspace engine, and the subspace and
%   full-space engines agree (single orientation, ideal pulses, echo
%   pathways, at most 64 points per kind-1 dimension, i.e. no binning)

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
kinds{1} = 3;

% three incremented intervals in one dimension
Exps{2} = Exp;
Exps{2}.Sequence = {p90 0.12 p90 0.1 p90 0.12};
Exps{2}.nPoints = 80;
Exps{2}.Dim1 = {'d1,d2,d3' 0.01};
Exps{2}.Pathways = '+0-';
kinds{2} = 3;

% two dimensions incrementing the same intervals, each a binning scheme on
% its own: the second one (more points) stays kind 1
Exps{3} = Exp;
Exps{3}.Sequence = {p90 0.12 p90 0.1 p90 0.12};
Exps{3}.nPoints = [6 40];
Exps{3}.Dim1 = {'d1,d3' 0.04};
Exps{3}.Dim2 = {'d1,d3' 0.01};
Exps{3}.Pathways = '+0-';
kinds{3} = [3 1];

tol = 1e-5;  % high-field approximation of the subspace engine
for k = 1:numel(Exps)
  E = Exps{k};
  [seq,reasons] = runprivate('sf_parsesequence',E,false);
  ok(k,1) = isequal(seq.kind,kinds{k}) && isempty(reasons); %#ok<AGROW>
  [~,ySub,info] = saffron(Sys,E);
  ySubs{k} = ySub; %#ok<AGROW>
  ok(k,2) = strcmp(info.Method,'subspace');
  E.DetWindow = 0;
  [~,yFull] = saffron(Sys,E,struct('Method','fullspace'));
  ok(k,3) = areequal(ySub,yFull,tol,'rel');
  % the signal is modulated much more than the tolerance
  ok(k,4) = std(real(ySub(:)))/max(abs(ySub(:))) > 10*tol;
end
ok = ok(:);

% Split delays, point by point, vs 2p ESEEM with the summed delays
E = Exp;
E.Sequence = {p90 0.15 p180 0.15};
E.nPoints = 64;
E.Dim1 = {'d1,d2' 0.02};
E.Pathways = '+-';
[~,y2p] = saffron(Sys,E);
ok(end+1) = areequal(ySubs{1},y2p,1e-12,'rel');

end
