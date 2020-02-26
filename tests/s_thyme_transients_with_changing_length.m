function ok = test()

% this test asserts the output is correctly returned for the various cases
% where the length of the transient (in points) changes due to events
% changing their length.
Sys.ZeemanFreq = 9.400;

Pulse.tp = 0.030;
Pulse.Flip = pi;
Pulse.Frequency =  200;

% Experiment/Sequence
tau = 0.3;
Exp.Sequence = {tau Pulse tau}; % us
Exp.mwFreq = 9.400; % GHz

% Options

Exp.Dim1 = {'p1.Position' 0.1};
Exp.Dim2 = {'p1.Position' 0.01};

Opt.IntTimeStep = 3.8e-6;

TimeAxisChanged = zeros(1,4);
SignalChanged = zeros(1,4);

% 1st test, 1 indirect dimension with one detection operator
Exp1 = Exp;
Exp1.nPoints = [3];
Exp1.DetOperator = {'z1'}; 

[TimeAxis1, Signal1] = spidyan(Sys,Exp1,Opt);

TimeAxisChanged(1) = any([~iscell(TimeAxis1), size(TimeAxis1) ~= [1 3], size(TimeAxis1{1}) ~= [1,165790], size(TimeAxis1{3}) ~= [1,165789]]);
SignalChanged(1) = any([~iscell(Signal1), size(Signal1) ~= [1 3], size(Signal1{1}) ~= [165790,1], size(Signal1{3}) ~= [165789,1]]);

% 2nd test, 1 indirect dimension with two detection operators
Exp2 = Exp;
Exp2.nPoints = [3];
Exp2.DetOperator = {'z1','z1'}; 

[TimeAxis2, Signal2] = spidyan(Sys,Exp2,Opt);

TimeAxisChanged(2) = any([~iscell(TimeAxis2), size(TimeAxis2) ~= [1 3], size(TimeAxis2{1}) ~= [1,165790], size(TimeAxis2{3}) ~= [1,165789]]);
SignalChanged(2) = any([~iscell(Signal2), size(Signal2) ~= [1 3], size(Signal2{1}) ~= [165790,2], size(Signal2{3}) ~= [165789,2]]);

% 3rd test, 2 indirect dimensions with one detection operator
Exp3 = Exp;
Exp3.nPoints = [3 2];
Exp3.DetOperator = {'z1'}; 

[TimeAxis3, Signal3] = spidyan(Sys,Exp3,Opt);

TimeAxisChanged(3) = any([~iscell(TimeAxis3), size(TimeAxis3) ~= [3 2], size(TimeAxis3{1}) ~= [1,165790], size(TimeAxis3{3}) ~= [1,165789]]);
SignalChanged(3) = any([~iscell(Signal3), size(Signal3) ~= [3 2], size(Signal3{1}) ~= [165790,1], size(Signal3{3}) ~= [165789,1]]);

% 4th test, 2 indirect dimensions with two detection operators
Exp4 = Exp;
Exp4.nPoints = [3 2];
Exp4.DetOperator = {'z1' 'z1'}; 

[TimeAxis4, Signal4] = spidyan(Sys,Exp4,Opt);

TimeAxisChanged(4) = any([~iscell(TimeAxis4), size(TimeAxis4) ~= [3 2], size(TimeAxis4{1}) ~= [1,165790], size(TimeAxis4{3}) ~= [1,165789]]);
SignalChanged(4) = any([~iscell(Signal4), size(Signal4) ~= [3 2], size(Signal4{1}) ~= [165790,2], size(Signal4{3}) ~= [165789,2]]);

ok(1) = all(~TimeAxisChanged);
ok(2) = all(~SignalChanged);
