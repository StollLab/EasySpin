function [err,data] = test(opt,olddata)
% Check pulse() resonator simulation and compensation
%--------------------------------------------------------------------------

% Rectangular pulse at resonant frequency
%--------------------------------------------------------------------------
Par.tp = 0.200; % us
[tpulse,IQ] = pulse(Par);

% Transfer function
Par.mwFreq = 9.50; % GHz
Par.ResonatorFrequency = 9.50; % GHz
Par.ResonatorQL = 500;

Opt.Resonator = 'compensate';
[tcomp,IQcomp] = pulse(Par,Opt);

% Test pulse function
Par_.tp = tcomp(end);
Par_.IQ = IQcomp;

Opt.Resonator = 'simulate';
Par_.ResonatorFrequency = Par.ResonatorFrequency;
Par_.ResonatorQL = Par.ResonatorQL;
Par_.mwFreq = Par.mwFreq;

[ttest,IQtest] = pulse(Par_,Opt);

% Reposition
splineinterp = @(x) interp1(ttest+x,IQtest,tpulse,'spline',0);
fitfunc = @(x) sqrt(sum((splineinterp(x)-IQ).^2)/numel(IQ));
x0 = 0;
x = fminsearch(fitfunc,x0);
IQtest = interp1(ttest+x,IQtest,tpulse,'spline',0);

suberr(1) = ~areequal(IQtest(2:end-1),IQ(2:end-1),0.01*max(IQ));


% Rectangular pulse shifted from resonator center
%--------------------------------------------------------------------------
clearvars -except suberr
Par.tp = 0.200; % us
[tpulse,IQ] = pulse(Par);

% Transfer function
Par.mwFreq = 15.00; % GHz
Par.ResonatorFrequency = 15.50; % GHz
Par.ResonatorQL = 500;

Opt.Resonator = 'compensate';
Par.TimeStep = 0.00025;
[tcomp,IQcomp] = pulse(Par,Opt);

% Test pulse function
Par_.tp = tcomp(end);
Par_.IQ = IQcomp;

Opt.Resonator = 'simulate';
Par_.ResonatorFrequency = Par.ResonatorFrequency;
Par_.ResonatorQL = Par.ResonatorQL;
Par_.mwFreq = Par.mwFreq;

[ttest,IQtest] = pulse(Par_,Opt);

% Reposition
splineinterp = @(x) interp1(ttest+x,IQtest,tpulse,'spline',0);
fitfunc = @(x) sqrt(sum((splineinterp(x)-IQ).^2)/numel(IQ));
x0 = 0;
x = fminsearch(fitfunc,x0);
IQtest = interp1(ttest+x,IQtest,tpulse,'spline',0);

suberr(2) = ~areequal(IQtest(2:end-1),IQ(2:end-1),0.01*max(IQ));


% Gaussian pulse off resonance
%--------------------------------------------------------------------------
clearvars -except suberr
Par.tp = 0.200; % us
Par.Type = 'gaussian';
Par.tFWHM = 0.050; % us
[tpulse,IQ] = pulse(Par);

% Transfer function
Par.mwFreq = 33.750; % GHz
Par.ResonatorFrequency = 33.800; % GHz
Par.ResonatorQL = 800;

Opt.Resonator = 'compensate';
[tcomp,IQcomp] = pulse(Par,Opt);

% Test pulse function
Par_.tp = tcomp(end);
Par_.IQ = IQcomp;

Opt.Resonator = 'simulate';
Par_.ResonatorFrequency = Par.ResonatorFrequency;
Par_.ResonatorQL = Par.ResonatorQL;
Par_.mwFreq = Par.mwFreq;

[ttest,IQtest] = pulse(Par_,Opt);

% Reposition
splineinterp = @(x) interp1(ttest+x,IQtest,tpulse,'spline',0);
fitfunc = @(x) sqrt(sum((splineinterp(x)-IQ).^2)/numel(IQ));
x0 = 0;
x = fminsearch(fitfunc,x0);
IQtest = interp1(ttest+x,IQtest,tpulse,'spline',0);

suberr(3) = ~areequal(IQtest,IQ,0.01*max(IQ));


% Sech/tanh pulse shifted from resonator frequency
%--------------------------------------------------------------------------
clearvars -except suberr

Par.tp = 0.200; % us
Par.Type = 'sech/tanh';
Par.beta = 10;
Par.Frequency = [-50 50]; % MHz
[tpulse,IQ] = pulse(Par);

% Transfer function
Par.mwFreq = 9.45; % GHz
Par.ResonatorFrequency = 9.50; % GHz
Par.ResonatorQL = 500;

Opt.Resonator = 'compensate';
[tcomp,IQcomp] = pulse(Par,Opt);

% Test pulse function
Par_.tp = tcomp(end);
Par_.IQ = IQcomp;
Par_.TimeStep = tpulse(2)-tpulse(1);

Opt.Resonator = 'simulate';
Par_.ResonatorFrequency = Par.ResonatorFrequency;
Par_.ResonatorQL = Par.ResonatorQL;
Par_.mwFreq = Par.mwFreq;

[ttest,IQtest] = pulse(Par_,Opt);

% Reposition
splineinterp = @(x) interp1(ttest+x,IQtest,tpulse,'spline');
fitfunc = @(x) sqrt(sum((splineinterp(x)-IQ).^2)/numel(IQ));
x0 = 0;
x = fminsearch(fitfunc,x0);
IQtest = interp1(ttest+x,IQtest,tpulse,'spline');

suberr(4) = ~areequal(IQtest,IQ,0.01*max(IQ));


err = any(suberr);

data = [];

