function ok = test() %#ok<FNDEF>

% Subspace engine with finite pulses: Opt.nOffsets = 1 gives a single offset
% of zero, the same as a vanishingly narrow inhomogeneous line
% (Sys.lwAnisoG)

Sys.Nucs = '1H';
Sys.A = [2 6];  % MHz

Exp.Field = 340;  % mT
p90 = struct('Flip',pi/2,'tp',0.004);
p180 = struct('Flip',pi,'tp',0.008);
Exp.Sequence = {p90 0.15 p180 0.15};
Exp.Dim1 = {'d1,d2' 0.02};
Exp.nPoints = 4;

Opt.nOffsets = 1;
[~,y1] = saffron(Sys,Exp,Opt);

Opt.nOffsets = 3;
Sys.lwAnisoG = 1e-9*[1 1 1];  % MHz
[~,y3] = saffron(Sys,Exp,Opt);

ok = areequal(y1,y3,1e-10,'rel');
