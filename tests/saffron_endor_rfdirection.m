function ok = test() %#ok<FNDEF>

% Direction of the RF field. The RF field is linearly polarized along
% y(Lab). The powder grid has orientations (phi,theta) with chi = 0, and
% rotating a molecule about the field by chi is equivalent to rotating the
% RF field by -chi. So, in powders, both engines average the signal over
% the RF field direction in the plane perpendicular to the field.
% - Powder (axial system, Dinfh grid), with and without the hyperfine
%   enhancement: equals the weighted sum of crystal orientations
%   (0,theta_k), each averaged explicitly over the rotation about the field
%   (Exp.SampleRotation).
% - Without enhancement near the cancellation condition A = 2*nuI, the RF
%   direction matters strongly: the RF field along y(Lab) only would be
%   far off.
% - Crystals: the RF field is along y(Lab), so the signal depends on the
%   rotation about the field, in the same way in both engines.
% - Strong RF pulses at low frequencies (full engine): the rotating-wave
%   error depends on the RF field direction, and RF pulses are stepped
%   where it exceeds 1%. The stepped points are the same for all
%   directions, so the signal is smooth in the direction and the average
%   converges with few directions (with a choice per direction, it did
%   not converge with 256). It agrees with the lab frame, where all RF
%   pulses are stepped, within the error of the unstepped points (about
%   2% of the largest ENDOR effect).

p90.Flip = pi/2;
Sys.Nucs = '1H';
Sys.A = [-2 -2 8];  % MHz, axial: the RF direction matters off the axes
rf.Channel = 'rf';
rf.tp = 2;  % µs
rf.Flip = pi;
Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.Sequence = {p90 0.2 p90 0 rf 0.1 p90 0.2};
Exp.Pathways = '+00-';
Exp.DetPhase = -pi/2;
Exp.rfFreq = 10;  % MHz
Exp.nPoints = 51;
Exp.Dim1 = {'rfFreq',0.2};

Opt.GridSymmetry = 'Dinfh';
Opt.GridFrame = [0 0 0];
Opt.GridSize = 7;
grid = sphgrid(Opt.GridSymmetry,Opt.GridSize);
nRot = 32;
rho = (0:nRot-1)*pi/nRot;  % rotations about the field, half a turn (period pi)

methods = {'manifold','full'};
ExpC = Exp;
ExpC.CrystalSymmetry = 1;
ok = false(1,0);
for m = 1:2
  if m==2
    % full engine: a few RF frequencies around one ENDOR line
    Exp.rfFreq = 17.4;  % MHz
    Exp.nPoints = 5;
    ExpC.rfFreq = Exp.rfFreq;
    ExpC.nPoints = Exp.nPoints;
    nRot = 16;
    rho = (0:nRot-1)*pi/nRot;
  end
  for enh = [true false]
    OptM = Opt;
    OptM.Method = methods{m};
    OptM.rfEnhancement = enh;
    [~,yPowder,info] = saffron(Sys,Exp,OptM);
    yCrystals = rotationAverage(Sys,ExpC,grid,rho,struct('Method',methods{m},'rfEnhancement',enh));
    ok(end+1) = areequal(yPowder,yCrystals,1e-5*max(abs(yPowder)),'abs') && ...
      info.nRFDirections>=4; %#ok<AGROW>
  end
end

% Without hyperfine enhancement, near the cancellation condition A = 2*nuI
% (manifold engine): the tilt of the nuclear quantization axes makes the
% RF direction matter strongly. The powder equals the rotation average, and
% differs clearly from the crystal sum with the RF field along y(Lab) only.
SysC = Sys;
SysC.A = [24 24 36];  % MHz, nuI = 14.9 MHz
ExpN = Exp;
ExpN.rfFreq = 0.5;  % MHz
ExpN.nPoints = 41;
ExpN.Dim1 = {'rfFreq',1};
ExpNC = ExpN;
ExpNC.CrystalSymmetry = 1;
OptN = struct('Method','manifold','rfEnhancement',false);
OptP = Opt;
OptP.Method = 'manifold';
OptP.rfEnhancement = false;
[~,yPowder] = saffron(SysC,ExpN,OptP);
yAverage = rotationAverage(SysC,ExpNC,grid,(0:31)*pi/32,OptN);
yAlongY = rotationAverage(SysC,ExpNC,grid,0,OptN);
ExpOff = ExpNC;
ExpOff.Sequence{5}.Flip = 0;
ExpOff.nPoints = 1;
yOff = rotationAverage(SysC,ExpOff,grid,0,OptN);
endorMax = max(abs(yOff-yAverage));
ok(end+1) = areequal(yPowder,yAverage,1e-3*endorMax,'abs') && ...
  max(abs(yAlongY-yAverage))>0.05*endorMax;

% Crystals: dependence on the rotation about the field, in both engines
ExpC = Exp;
ExpC.CrystalSymmetry = 1;
ExpC.MolFrame = [0 pi/4 0];
dy = cell(1,2);
for m = 1:2
  ExpC.SampleRotation = {[0 0 1],0};
  [~,y0] = saffron(Sys,ExpC,struct('Method',methods{m}));
  ExpC.SampleRotation = {[0 0 1],pi/2};
  [~,y90] = saffron(Sys,ExpC,struct('Method',methods{m}));
  dy{m} = y90 - y0;
end
ok(end+1) = max(abs(dy{2}))>1e-3*max(abs(y0)) && ...
  areequal(dy{1},dy{2},0.05*max(abs(dy{2})),'abs');

% Strong RF pulses at low frequencies, without enhancement (fast in the
% lab frame)
SysS.Nucs = '1H';
SysS.A = [20 20 36];  % MHz, lines down to 2 MHz
ExpS = Exp;
ExpS.Sequence{4} = 0.1;  % µs
ExpS.Sequence{5} = struct('Channel','rf','tp',2,'B2',2);  % mT
ExpS.rfFreq = 2;  % MHz
ExpS.nPoints = 11;
ExpS.Dim1 = {'rfFreq',0.5};
ExpOff = ExpS;
ExpOff.Sequence{5}.B2 = 0;
ExpOff.nPoints = 1;
OptS = Opt;
OptS.Method = 'full';
OptS.rfEnhancement = false;
OptL = OptS;
OptL.SimFreq = 0;
[~,yS,infoS] = saffron(SysS,ExpS,OptS);
[~,ySoff] = saffron(SysS,ExpOff,OptS);
[~,yL] = saffron(SysS,ExpS,OptL);
[~,yLoff] = saffron(SysS,ExpOff,OptL);
endorL = real(yLoff-yL);
ok(end+1) = infoS.nRFDirections<=16 && ...
  areequal(real(ySoff-yS),endorL,0.03*max(abs(endorL)),'abs');

end

% Weighted sum of crystal orientations (0,theta_k) of the powder grid, each
% averaged over the rotations rho about the field (Exp.SampleRotation)
function y = rotationAverage(Sys,ExpC,grid,rho,Opt)
y = 0;
for k = 1:numel(grid.theta)
  ExpC.MolFrame = [0 grid.theta(k) 0];
  ExpC.SampleRotation = {[0 0 1],rho};
  [~,yk] = saffron(Sys,ExpC,Opt);
  y = y + grid.weights(k)*yk;
end
y = y/sum(grid.weights);
end
