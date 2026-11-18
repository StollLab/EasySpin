function ok = test()  %#ok<FNDEF>

% info.Pathways and info.EchoTimes: the contributing coherence pathways and
% the times of their echoes after the end of the sequence, sorted by echo time
% - all pathways: echo times, NaN for pathways without echo (listed last)
% - phase cycle: cancelled pathways are not listed, and the listed pathways
%   as Exp.Pathways give the same signal
% - 'echo' with a transient: the echoes within the window (full-space)
% - predefined experiment: its pathways
% - full-space engine: the same list for S = 1/2, empty for S = 1
% - finite pulses: user-defined waveforms, pulses switched on and off

Sys.g = 2;
p.Flip = pi/2;

Exp.mwFreq = 9.5;  % GHz
Exp.Field = mhz2mt(9500,2);  % mT
Exp.MolFrame = [0 0 0];
Exp.CrystalSymmetry = 1;
Exp.Sequence = {p 0.1 p 0.5 p};

% All pathways (tau = 0.1, T = 0.5): echo times from the delays, NaN without
% +1 coherence
[~,~,info] = saffron(Sys,Exp);
echoRows = ['+--'; '+0-'; '-+-'; '0+-'; '++-'];
noEchoRows = ['00-'; '0--'; '-0-'; '---'];
ok(1) = isequal(info.Pathways(1:5,:),echoRows) && ...
  areequal(info.EchoTimes(1:5),[-0.4; 0.1; 0.4; 0.5; 0.6],1e-12,'abs');
ok(2) = isequal(sortrows(info.Pathways(6:end,:)),sortrows(noEchoRows)) && ...
  all(isnan(info.EchoTimes(6:end)));

% Phase cycle on the first pulse: removes the pathways without coherence
% after it; the remaining pathways give the same signal as the phase cycle
Ec = Exp;
Ec.PhaseCycle = [0 0 0 0; pi 0 0 pi];
[~,yCycle,info] = saffron(Sys,Ec);
ok(3) = isequal(sortrows(info.Pathways),sortrows(['+0-'; '++-'; '+--'; '-+-'; '-0-'; '---']));
Ee = Exp;
Ee.Pathways = info.Pathways;
[~,yExplicit] = saffron(Sys,Ee);
ok(4) = abs(yCycle)>0.1 && areequal(yExplicit,yCycle,1e-10,'rel');

% Echoes within a transient window (full-space engine)
Et = Exp;
Et.DetWindow = [0 2];  % µs
Et.DetTimeStep = 0.01;  % µs
Et.Pathways = 'echo';
[~,~,info] = saffron(Sys,Et);
ok(5) = strcmp(info.Method,'fullspace') && ...
  isequal(info.Pathways,['+0-'; '-+-'; '0+-'; '++-']) && ...
  areequal(info.EchoTimes,[0.1; 0.4; 0.5; 0.6],1e-12,'abs');

% Predefined experiment: its pathway
Eh = rmfield(Exp,'Sequence');
Eh.Sequence = 'HYSCORE';
Eh.tau = 0.1;  % µs
Eh.dt = 0.02;  % µs
Eh.nPoints = [4 4];
[~,~,info] = saffron(Sys,Eh);
ok(6) = isequal(info.Pathways,'+00-') && areequal(info.EchoTimes,0,1e-12,'abs');

% Full-space engine: same list as the subspace engine for S = 1/2
[~,~,infoSub] = saffron(Sys,Exp);
[~,~,infoFull] = saffron(Sys,Exp,struct('Method','fullspace'));
ok(7) = isequal(infoFull.Pathways,infoSub.Pathways) && ...
  isequaln(infoFull.EchoTimes,infoSub.EchoTimes);

% Finite pulses: user-defined waveforms and rectangular pulses give the same
% list (no pruning by transfer amplitudes)
t = 0:0.001:0.016;  % µs
pIQ = struct('IQ',ones(size(t))/(4*0.016),'t',t);  % pi/2 pulse, amplitude in MHz
pR = struct('Flip',pi/2,'tp',0.016);  % µs
Ef = Exp;
Ef.Sequence = {pIQ 0.1 pIQ 0.1};
Opt = struct('IntTimeStep',1e-4);  % µs
[~,~,infoIQ] = saffron(Sys,Ef,Opt);
Ef.Sequence = {pR 0.1 pR 0.1};
[~,~,infoRect] = saffron(Sys,Ef,setfield(Opt,'Method','fullspace')); %#ok<SFLD>
ok(end+1) = isequal(infoIQ.Pathways,['+-'; '0-'; '--']) && isequal(infoRect.Pathways,infoIQ.Pathways);

% A finite pulse switched on and off: the pathways that the phase cycle keeps
% (no phase cycle: all 9 combinations of +, 0, - after the first two pulses)
Eo = Exp;
Eo.Sequence = {pR 0.1 struct('Flip',pi,'tp',0.032) 0.2 pR 0.1};
Eo.OnOff = 2;
[~,y,info] = saffron(Sys,Eo,Opt);
ok(end+1) = strcmp(info.Method,'fullspace') && abs(y)>0.1 && size(info.Pathways,1)==9;

% Full-space engine with S = 1: no pathway list
Sys1.S = 1;
Sys1.g = 2;
pA = struct('tp',0.01,'Amplitude',25);  % µs, MHz
E1 = Exp;
E1.Sequence = {pA 0.1 pA 0.5 pA};
[~,~,info] = saffron(Sys1,E1);
ok(end+1) = strcmp(info.Method,'fullspace') && isempty(info.Pathways) && isempty(info.EchoTimes);

end
