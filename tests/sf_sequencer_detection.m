function ok = test()

Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % µs
Pulse.Flip = pi;
Pulse.Frequency = 1000* [-0.100 0.100];
Pulse.tp = 0.1;

Pulse2 = Pulse;
Pulse2.tp = 0.2;

Exp.Sequence = {Pulse 0.5 Pulse Pulse2};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % µs

Exp.mwFreq = 33.5;
Exp.DetSequence = [1 0 0 0]; 


% First Test ---------------------------------------

Opt.DetEvents = true;

[Events1] = runprivate('sf_buildevents',Exp,Opt);

Opt.DetEvents = [1 1 1 1];

[Events2] = runprivate('sf_buildevents',Exp,Opt);

% Second Test ---------------------------------------

Opt.DetEvents = false;

[Events3] = runprivate('sf_buildevents',Exp,Opt);

Opt.DetEvents = [0 0 0 0];

[Events4] = runprivate('sf_buildevents',Exp,Opt);

Opt.DetEvents = [];

[Events5] = runprivate('sf_buildevents',Exp,Opt);

Opt = rmfield(Opt,'DetEvents');

[Events6] = runprivate('sf_buildevents',Exp,Opt);

% Third Test ---------------------------------------

Opt.DetEvents = [1 1];

[Events7] = runprivate('sf_buildevents',Exp,Opt);

Opt.DetEvents = [1 1 0 0];

[Events8] = runprivate('sf_buildevents',Exp,Opt);


ok = all([isequal(Events1,Events2) isequal(Events3,Events4) ...
     isequal(Events3,Events4) isequal(Events3,Events5) ...
     isequal(Events3,Events6) isequal(Events7,Events8)]);
