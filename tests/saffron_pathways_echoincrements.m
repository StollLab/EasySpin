function ok = test()  %#ok<FNDEF>

% Exp.Pathways = 'echo' selects the pathways whose echo forms within the
% detection window (Exp.DetWindow) at all points of the experiment, for any
% increments and decrements of the delays
% - equal steps, unequal steps, decrements, two dimensions, values that do
%   not start at 0, non-linear values, several incremented delays in one
%   interval, moving pulses
% - echo position in Exp.DetWindow, transient windows, echoes that move
%   within the window, finite pulses (window widened by the pulse lengths,
%   also for pulses that are finite only at some points)
% - both engines use the same pathways at all points
% - errors: no echo within the window at all points, DetWindow = 'sequence',
%   increments of delays before the first pulse (subspace engine only)

p.Flip = pi/2;

% Stimulated echo with equal steps
E = struct('Sequence',{{p 0.2 p 0.2 p 0.2}},'nPoints',4);
E.Dim1 = {'d1,d3' 0.02};
ok(1) = isequal(echopathways(E),[1 0 -1]);

% Unequal steps: [0 1 -1] refocuses only at the first point
Eu = E;
Eu.Dim1 = {'d1' 0.02; 'd2' 0.04; 'd3' 0.02};
ok(2) = isequal(echopathways(Eu),[1 0 -1]);

% Decrement: echo position stays constant
Ed = struct('Sequence',{{p 0.1 p 0.1 p 0.2}},'nPoints',4);
Ed.Dim1 = {'d1' 0.02; 'd2' -0.02};
ok(3) = isequal(echopathways(Ed),[1 1 -1]);

% Two dimensions
E2 = struct('Sequence',{{p 0.2 p 0.2 p 0.2}},'nPoints',[4 3]);
E2.Dim1 = {'d1,d3' 0.02};
E2.Dim2 = {'d2' 0.02};
ok(4) = isequal(echopathways(E2),[1 0 -1]);

% Values that do not start at 0
Eo = struct('Sequence',{{p 0.1 p 0.2}},'nPoints',3);
Eo.Dim1 = {'d1' [0.1 0.12 0.14]; 'd2' [0 0.02 0.04]};
ok(5) = isequal(echopathways(Eo),[1 -1]);

% Several incremented delays in one interval each count
Es = struct('Sequence',{{p 0.1 p 0.1 p 0.1 0.1}},'nPoints',4);
Es.Dim1 = {'d1,d2,d3,d4' 0.02};
ok(6) = isequal(echopathways(Es),[1 1 -1]);

% Moving pulse: the delays before and after it change in opposite directions
Em = struct('Sequence',{{p 0.1 p 0.1 p 0.2}},'nPoints',4);
Em.Dim1 = {'p2.Position' 0.02};
ok(7) = isequal(echopathways(Em),[1 1 -1]);

% Incremented RF pulse length, with both names of the field
r = struct('Channel','rf','tp',0.02,'Flip',pi);
Er = struct('Sequence',{{p 0.1 r 0.1 p 0.22}},'nPoints',4,'rfFreq',14);
for field = {'tp','t'}
  Er.Dim1 = {['p2.' field{1}] 0.01; 'd3' 0.01};
  ok(end+1) = isequal(echopathways(Er),[1 1 -1]); %#ok<AGROW>
end

% A delay before the first pulse has no effect
Eb = struct('Sequence',{{0.3 p 0.2 p 0.2}},'nPoints',4);
Eb.Dim1 = {'d1' 0.05; 'd2' 0.02; 'd3' 0.02};
ok(end+1) = isequal(echopathways(Eb),[1 -1]);

% Echo moves away from the detection point
noEcho = {};
noEcho{1} = struct('Sequence',{{p 0.2 p 0.2}},'nPoints',4);
noEcho{1}.Dim1 = {'d1' 0.02; 'd2' -0.02};
noEcho{2} = struct('Sequence',{{p 0.1 0.1 p 0.2}},'nPoints',4);
noEcho{2}.Dim1 = {'d1,d2,d3' 0.02};
noEcho{3} = struct('Sequence',{{p 0.2 p 0.2}},'nPoints',4);
noEcho{3}.Dim1 = {'p2.Position' 0.02};
for k = 1:numel(noEcho)
  ok(end+1) = throwsError(noEcho{k},'No pathway has its echo');  %#ok<AGROW>
end

% Non-linear delay values
En = struct('Sequence',{{p 0.2 p 0.2}},'nPoints',4);
En.Dim1 = {'d1,d2' [0 0.05 0.2 0.21]};
ok(end+1) = isequal(echopathways(En),[1 -1]);

% Echo position in the detection window: the sequence ends with a pulse
Ep = struct('Sequence',{{p 0.1 p 0.5 p}},'DetWindow',0.1);
ok(end+1) = isequal(echopathways(Ep),[1 0 -1]);

% Transient: all echoes within the window, at tau (+0-), T-tau (-+-), T (0+-)
% and T+tau (++-) after the last pulse; +-- forms its echo before the last pulse
Et = Ep;
Et.DetWindow = [0 2];
ok(end+1) = isequal(echopathways(Et),sortrows([1 0 -1; -1 1 -1; 0 1 -1; 1 1 -1]));
Et.DetWindow = [0 0.45];
ok(end+1) = isequal(echopathways(Et),sortrows([1 0 -1; -1 1 -1]));

% Echo that moves within the window, and out of it
Ew = struct('Sequence',{{p 0.2 p 0.2}},'nPoints',4,'DetWindow',[0 0.1]);
Ew.Dim1 = {'d1' 0.02};  % echo at 0, 0.02, 0.04, 0.06 after the end
ok(end+1) = isequal(echopathways(Ew),[1 -1]);
Ew.DetWindow = [0 0.05];
ok(end+1) = throwsError(Ew,'No pathway has its echo');

% Finite pulses: the echo maximum is shifted by about 2*tp/pi, within the
% total pulse length; with ideal pulses, the shifted point is not an echo
pf = struct('Flip',pi/2,'tp',0.016);
Ef = struct('Sequence',{{pf 0.2 pf}},'DetWindow',0.2+2*pf.tp/pi);
ok(end+1) = isequal(echopathways(Ef),[1 -1]);
Ef.Sequence = {p 0.2 p};
ok(end+1) = throwsError(Ef,'No pathway has its echo');

% A pulse that is ideal at the first point and finite later widens the
% window by its largest length
Ef0 = struct('Sequence',{{p 0.2 p 0.2}},'nPoints',3,'DetWindow',0.015);
Ef0.Dim1 = {'p1.tp' 0.01};  % 0, 10 and 20 ns
ok(end+1) = isequal(echopathways(Ef0),[1 -1]);

% Detection during the whole sequence: no detection window
Eq = struct('Sequence',{{p 0.2 p 0.2}},'DetWindow','sequence');
ok(end+1) = throwsError(Eq,'DetWindow = ''sequence''');

% Both engines: same signal as with the explicit pathway
Sys.g = [2.00 2.01 2.03];
Exp = Eu;
Exp.mwFreq = 9.5;  % GHz
Exp.Field = 338;  % mT
Exp.MolFrame = [0.2 0.8 0.1];
Exp.CrystalSymmetry = 1;
Exp.DetWindow = 0;
tol = 1e-10;
for method = {'subspace','fullspace'}
  Opt.Method = method{1};
  Exp.Pathways = 'echo';
  [~,yEcho] = saffron(Sys,Exp,Opt);
  Exp.Pathways = '+0-';
  [~,yExplicit] = saffron(Sys,Exp,Opt);
  ok(end+1) = areequal(yEcho,yExplicit,tol,'rel'); %#ok<AGROW>
end

% Echo position in the detection window, in both engines (the subspace
% engine by default)
Exp = Ep;
Exp.mwFreq = 9.5;  % GHz
Exp.Field = 338;  % mT
Exp.MolFrame = [0.2 0.8 0.1];
Exp.CrystalSymmetry = 1;
for method = {'auto','fullspace'}
  Opt.Method = method{1};
  Exp.Pathways = 'echo';
  [~,yEcho,info] = saffron(Sys,Exp,Opt);
  Exp.Pathways = '+0-';
  [~,yExplicit] = saffron(Sys,Exp,Opt);
  ok(end+1) = strcmp(info.Method,strrep(method{1},'auto','subspace')) && abs(yEcho)>0.1 && ...
    areequal(yEcho,yExplicit,tol,'rel'); %#ok<AGROW>
end

% Incrementing only a delay before the first pulse: no effect on the initial
% state of the subspace engine (error); allowed in the full-space engine
Exp = rmfield(Eb,'Dim1');
Exp.Dim1 = {'d1' 0.05};
Exp.mwFreq = 9.5;  % GHz
Exp.Field = 338;  % mT
Exp.MolFrame = [0.2 0.8 0.1];
Exp.CrystalSymmetry = 1;
try
  saffron(Sys,Exp);
  ok(end+1) = false;
catch ME
  ok(end+1) = contains(ME.message,'before the first pulse has no effect');
end
[~,y] = saffron(Sys,Exp,struct('Method','fullspace'));
ok(end+1) = numel(y)==4 && areequal(y,y(1)*ones(size(y)),1e-10,'rel');

end

function co = echopathways(Exp)
Exp.Pathways = 'echo';
seq = runprivate('sf_parsesequence',Exp,false);
co = unique(seq.co,'rows');  % alpha and beta populations both have order 0
end

function thrown = throwsError(Exp,text)
thrown = false;
try
  echopathways(Exp);
catch ME
  thrown = contains(ME.message,text);
end
end
