function ok = test() %#ok<FNDEF>

% Exp.PhaseCycle
% - the order of the steps does not matter
% - full-space method: phases are compared modulo 2*pi and with a
%   tolerance
% - errors for wrong input

Sys.Nucs = '1H';
Sys.A = [2 3 9];  % MHz

p90.Flip = pi/2;  % rad

Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.MolFrame = [0 0.4 0.5];
Exp.CrystalSymmetry = 1;
Exp.Sequence = {p90 0.12 p90 0.3 p90 0.12};
Exp.DetPhase = -pi/2;  % rad

[ph1,ph3] = ndgrid((0:3)*pi/2,[0 pi]);
PC = [ph1(:) 0*ph1(:) ph3(:) -ph1(:)-ph3(:)];

tol = 1e-10;
methods = {'subspace','fullspace'};
for m = 1:2
  Opt.Method = methods{m};

  E = Exp;
  E.PhaseCycle = PC;
  [~,y] = saffron(Sys,E,Opt);

  % order of the steps
  E.PhaseCycle = PC(end:-1:1,:);
  [~,yRev] = saffron(Sys,E,Opt);
  ok(m) = areequal(yRev,y,tol,'rel'); %#ok<AGROW>
end

% Full-space method: phases that differ by 2*pi or by rounding errors are
% the same setting, so the cycle stays nested (a single propagation run)
p90f = struct('Flip',pi/2,'tp',0.01);
E = Exp;
E.Sequence = {p90f 0.12 p90f 0.3 p90f 0.12};
Opt = struct('Method','fullspace');
E.PhaseCycle = PC;
[~,y] = saffron(Sys,E,Opt);
PCmod = PC;
PCmod(PCmod(:,1)==3*pi/2,1) = -pi/2;
PCmod(PCmod(:,3)==pi,3) = 3*pi;
PCmod(:,1:3) = PCmod(:,1:3) + 1e-14*reshape(sin(1:24),8,3);
E.PhaseCycle = PCmod;
Opt.Verbosity = 1;
[logText,yMod] = evalc('saffron(Sys,E,Opt)');
ok(end+1) = areequal(yMod,y,1e-10,'rel') && contains(logText,'phase cycle: 1 propagation run(s)');

% Errors
errorCases = {
  {'PhaseCycle',[0 0; pi pi]}, '4 columns'
  {'PhaseCycle',{[0 1; pi -1]}}, 'no longer supported'
  };
for k = 1:size(errorCases,1)
  E = Exp;
  E.(errorCases{k,1}{1}) = errorCases{k,1}{2};
  try
    saffron(Sys,E);
    ok(end+1) = false; %#ok<AGROW>
  catch err
    ok(end+1) = contains(err.message,errorCases{k,2}); %#ok<AGROW>
  end
end

% Predefined experiments set their pathways themselves
E = struct('Field',350,'Sequence','MimsENDOR','tau',0.1,'rfFreq',[10 20],'PhaseCycle',PC);
try
  saffron(Sys,E);
  ok(end+1) = false;
catch err
  ok(end+1) = contains(err.message,'Exp.PhaseCycle cannot be used with predefined experiments');
end
