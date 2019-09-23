function [err,data] = test(opt,olddaa)
% Self-consistency check for resonator() 'simulate' and 'compensate'
% options
%--------------------------------------------------------------------------

% Rectangular pulse at resonant frequency
%--------------------------------------------------------------------------
Par.tp = 0.200; % us
[tpulse,IQ] = pulse(Par);

% Transfer function
mwFreq = 9.50; % GHz
ResonatorFrequency = 9.50; % GHz
ResonatorQL = 500;

% Resonator compensation
[tcomp,IQcomp] = resonator(tpulse,IQ,mwFreq,ResonatorFrequency,ResonatorQL,'compensate');

% Simulate effect of resonator on compensated pulse shape
[ttest,IQtest] = resonator(tcomp,IQcomp,mwFreq,ResonatorFrequency,ResonatorQL,'simulate');

% Reposition
splineinterp = @(x) interp1(ttest+x,IQtest,tpulse,'spline',0);
fitfunc = @(x) sqrt(sum((splineinterp(x)-IQ).^2)/numel(IQ));
x0 = 0;
x = fminsearch(fitfunc,x0);
IQtest = interp1(ttest+x,IQtest,tpulse,'spline',0);

suberr(1) = ~areequal(IQtest(2:end-1),IQ(2:end-1),0.01*max(IQ),'abs');


% Rectangular pulse shifted from resonator center
%--------------------------------------------------------------------------
clear Par Opt
Par.tp = 0.200; % us
[tpulse,IQ] = pulse(Par);

% Transfer function
mwFreq = 15.00; % GHz
ResonatorFrequency = 15.50; % GHz
ResonatorQL = 500;

% Resonator compensation
Opt.TimeStep = 0.00025;
[tcomp,IQcomp] = resonator(tpulse,IQ,mwFreq,ResonatorFrequency,ResonatorQL,'compensate',Opt);

% Simulate effect of resonator on compensated pulse shape
[ttest,IQtest] = resonator(tcomp,IQcomp,mwFreq,ResonatorFrequency,ResonatorQL,'simulate');

% Reposition
splineinterp = @(x) interp1(ttest+x,IQtest,tpulse,'spline',0);
fitfunc = @(x) sqrt(sum((splineinterp(x)-IQ).^2)/numel(IQ));
x0 = 0;
x = fminsearch(fitfunc,x0);
IQtest = interp1(ttest+x,IQtest,tpulse,'spline',0);

suberr(2) = ~areequal(IQtest(2:end-1),IQ(2:end-1),0.01*max(IQ),'abs');


% Gaussian pulse off resonance
%--------------------------------------------------------------------------
clear Par
Par.tp = 0.200; % us
Par.Type = 'gaussian';
Par.tFWHM = 0.050; % us
[tpulse,IQ] = pulse(Par);

% Transfer function
mwFreq = 33.750; % GHz
ResonatorFrequency = 33.800; % GHz
ResonatorQL = 800;

% Resonator compensation
[tcomp,IQcomp] = resonator(tpulse,IQ,mwFreq,ResonatorFrequency,ResonatorQL,'compensate');

% Simulate effect of resonator on compensated pulse shape
[ttest,IQtest] = resonator(tcomp,IQcomp,mwFreq,ResonatorFrequency,ResonatorQL,'simulate');

% Reposition
splineinterp = @(x) interp1(ttest+x,IQtest,tpulse,'spline',0);
fitfunc = @(x) sqrt(sum((splineinterp(x)-IQ).^2)/numel(IQ));
x0 = 0;
x = fminsearch(fitfunc,x0);
IQtest = interp1(ttest+x,IQtest,tpulse,'spline',0);

suberr(3) = ~areequal(IQtest,IQ,0.01*max(IQ),'abs');

% Sech/tanh pulse shifted from resonator frequency
%--------------------------------------------------------------------------
clear Par

Par.tp = 0.200; % us
Par.Type = 'sech/tanh';
Par.beta = 10;
Par.Frequency = [-50 50]; % MHz
[tpulse,IQ] = pulse(Par);

% Transfer function
mwFreq = 9.45; % GHz
ResonatorFrequency = 9.50; % GHz
ResonatorQL = 500;

% Resonator compensation
[tcomp,IQcomp] = resonator(tpulse,IQ,mwFreq,ResonatorFrequency,ResonatorQL,'compensate');

% Simulate effect of resonator on compensated pulse shape
[ttest,IQtest] = resonator(tcomp,IQcomp,mwFreq,ResonatorFrequency,ResonatorQL,'simulate');

% Reposition
splineinterp = @(x) interp1(ttest+x,IQtest,tpulse,'spline');
fitfunc = @(x) sqrt(sum((splineinterp(x)-IQ).^2)/numel(IQ));
x0 = 0;
x = fminsearch(fitfunc,x0);
IQtest = interp1(ttest+x,IQtest,tpulse,'spline');

suberr(4) = ~areequal(IQtest,IQ,0.01*max(IQ),'abs');

err = any(suberr);

data = [];

