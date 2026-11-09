function ok = test() %#ok<FNDEF>

% Exp.Pathways = 'echo' selects the pathways that refocus at all points of
% the experiment, for any increments and decrements of the delays
% - equal steps, unequal steps, decrements, two dimensions, values that do
%   not start at 0, non-linear values, several incremented delays in one
%   interval, moving pulses
% - both engines use the same pathways at all points
% - error: no pathway refocuses at all points

p.Flip = pi/2;

% Stimulated echo with equal steps
E = struct('Sequence',{{p 0.2 p 0.2 p 0.2}},'nPoints',4);
E.Dim1 = {'d1,d3' 0.02};
ok(1) = isequal(echoOrders(E),[1 0 -1]);

% Unequal steps: [0 1 -1] refocuses only at the first point
Eu = E;
Eu.Dim1 = {'d1' 0.02; 'd2' 0.04; 'd3' 0.02};
ok(2) = isequal(echoOrders(Eu),[1 0 -1]);

% Decrement: echo position stays constant
Ed = struct('Sequence',{{p 0.1 p 0.1 p 0.2}},'nPoints',4);
Ed.Dim1 = {'d1' 0.02; 'd2' -0.02};
ok(3) = isequal(echoOrders(Ed),[1 1 -1]);

% Two dimensions
E2 = struct('Sequence',{{p 0.2 p 0.2 p 0.2}},'nPoints',[4 3]);
E2.Dim1 = {'d1,d3' 0.02};
E2.Dim2 = {'d2' 0.02};
ok(4) = isequal(echoOrders(E2),[1 0 -1]);

% Values that do not start at 0
Eo = struct('Sequence',{{p 0.1 p 0.2}},'nPoints',3);
Eo.Dim1 = {'d1' [0.1 0.12 0.14]; 'd2' [0 0.02 0.04]};
ok(5) = isequal(echoOrders(Eo),[1 -1]);

% Several incremented delays in one interval each count
Es = struct('Sequence',{{p 0.1 p 0.1 p 0.1 0.1}},'nPoints',4);
Es.Dim1 = {'d1,d2,d3,d4' 0.02};
ok(6) = isequal(echoOrders(Es),[1 1 -1]);

% Moving pulse: the delays before and after it change in opposite directions
Em = struct('Sequence',{{p 0.1 p 0.1 p 0.2}},'nPoints',4);
Em.Dim1 = {'p2.Position' 0.02};
ok(7) = isequal(echoOrders(Em),[1 1 -1]);

% Incremented RF pulse length, with both names of the field
r = struct('Channel','rf','tp',0.02,'Flip',pi);
Er = struct('Sequence',{{p 0.1 r 0.1 p 0.22}},'nPoints',4,'rfFreq',14);
for field = {'tp','t'}
  Er.Dim1 = {['p2.' field{1}] 0.01; 'd3' 0.01};
  ok(end+1) = isequal(echoOrders(Er),[1 1 -1]); %#ok<AGROW>
end

% A delay before the first pulse has no effect
Eb = struct('Sequence',{{0.3 p 0.2 p 0.2}},'nPoints',4);
Eb.Dim1 = {'d1' 0.05; 'd2' 0.02; 'd3' 0.02};
ok(end+1) = isequal(echoOrders(Eb),[1 -1]);

% Echo moves away from the detection point
noEcho = {};
noEcho{1} = struct('Sequence',{{p 0.2 p 0.2}},'nPoints',4);
noEcho{1}.Dim1 = {'d1' 0.02; 'd2' -0.02};
noEcho{2} = struct('Sequence',{{p 0.1 0.1 p 0.2}},'nPoints',4);
noEcho{2}.Dim1 = {'d1,d2,d3' 0.02};
noEcho{3} = struct('Sequence',{{p 0.2 p 0.2}},'nPoints',4);
noEcho{3}.Dim1 = {'p2.Position' 0.02};
for k = 1:numel(noEcho)
  ok(end+1) = throwsError(noEcho{k},'No pathway refocuses'); %#ok<AGROW>
end

% Non-linear delay values
En = struct('Sequence',{{p 0.2 p 0.2}},'nPoints',4);
En.Dim1 = {'d1,d2' [0 0.05 0.2 0.21]};
ok(end+1) = isequal(echoOrders(En),[1 -1]);

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

end

function co = echoOrders(Exp)
Exp.Pathways = 'echo';
seq = runprivate('sf_parsesequence',Exp,struct);
co = unique(seq.co,'rows');  % alpha and beta populations both have order 0
end

function thrown = throwsError(Exp,text)
thrown = false;
try
  echoOrders(Exp);
catch ME
  thrown = contains(ME.message,text);
end
end
