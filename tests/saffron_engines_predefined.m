function ok = test() %#ok<FNDEF>

% The subspace and full-space engines give the same signals for the
% predefined experiments with ideal pulses (powder, S = 1/2, one 1H). The
% residual differences come from the frequency binning of the pathways
% engine (larger buffer expansion here) and the time step of the
% full-space engine.

Sys.Nucs = '1H';
Sys.A = [2 3 9];  % MHz
Sys.g = [2.0 2.01 2.02];

Base.Field = 350;  % mT
Base.mwFreq = 9.8;  % GHz
Base.tau = 0.12;  % µs
Base.dt = 0.02;  % µs

OptP.Method = 'subspace';
OptP.GridSize = 7;
OptP.Expand = 6;
OptT.Method = 'fullspace';
OptT.GridSize = 7;

names = {'2pESEEM','3pESEEM','4pESEEM','HYSCORE'};
tol = 5e-3;
for k = 1:numel(names)
  Exp = Base;
  Exp.Sequence = names{k};
  if k==4, Exp.nPoints = [8 8]; else, Exp.nPoints = 32; end
  [~,y1,info1] = saffron(Sys,Exp,OptP);
  [~,y2,info2] = saffron(Sys,Exp,OptT);
  ok(k) = areequal(y1,y2,tol,'rel') && strcmp(info1.Method,'subspace') && ...
    strcmp(info2.Method,'fullspace');
end
