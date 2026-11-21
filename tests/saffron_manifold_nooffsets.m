function ok = test() %#ok<FNDEF>

% Manifold engine with finite pulses: without inhomogeneous line, the
% signal is that of a single resonance offset of zero, the limit of a
% vanishingly narrow line (Sys.lwAnisoG)

Sys.Nucs = '1H';
Sys.A = [2 6];  % MHz

Exp.Field = 340;  % mT
p90 = struct('Flip',pi/2,'tp',0.004);
p180 = struct('Flip',pi,'tp',0.008);
Exp.Sequence = {p90 0.15 p180 0.15};
Exp.Dim1 = {'d1,d2' 0.02};
Exp.nPoints = 4;

[~,y0] = saffron(Sys,Exp);

Sys.lwAnisoG = 1e-9*[1 1 1];  % MHz
[~,yNarrow] = saffron(Sys,Exp);

ok = areequal(y0,yNarrow,1e-10,'rel');
