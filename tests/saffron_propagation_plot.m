function ok = test() %#ok<FNDEF>

% Plotting (saffron without output arguments) with one and with several
% detection operators, with and without an indirect dimension: one
% subplot per detection operator. A single point is printed, not plotted.

Sys.S = 1/2;

Pulse.Type = 'rectangular';
Pulse.tp = 0.02;  % µs
Pulse.Flip = pi/2;

Exp.Field = mhz2mt(9500);  % mT, resonance field for g = gfree at 9.5 GHz
Exp.mwFreq = 9.5;  % GHz
Exp.MolFrame = [0 0 0];  % single orientation
Exp.Sequence = {Pulse};
Exp.DetWindow = [0 0.02];  % µs

ops = {{'+1'},{'+1','x1','z1'}};
nAxes = zeros(2,2);
for iOp = 1:2
  for iDim = 1:2
    E = Exp;
    E.DetOperator = ops{iOp};
    if iDim==2
      E.nPoints = 3;
      E.Dim1 = {'p1.Flip' pi/4};
    end
    nAxes(iOp,iDim) = plotAxes(Sys,E);
  end
end
ok(1) = isequal(nAxes,[1 1; 3 3]);

% Single point: printed, no axes
E = Exp;
E.DetWindow = 0;
E.DetOperator = {'z1'};
[nAx,out] = plotAxes(Sys,E);
ok(2) = nAx==0 && contains(out,'signal');

end

% Number of axes in the figure from saffron without outputs (legends are
% not axes)
function [n,out] = plotAxes(Sys,Exp)
f = figure('Visible','off');
out = evalc('saffron(Sys,Exp)');
n = numel(findobj(f,'Type','axes'));
close(f);
end
