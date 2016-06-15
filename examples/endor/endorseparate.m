% 1H ENDOR, summed and separate output
%=====================================================================
% Demonstates the use of Options.Output
clear, clf

% Spin system
Spins.g = [2.25 2.25 2];
Spins = nucspinadd(Spins,'1H',[-1 2]*1 + 7);
Spins = nucspinadd(Spins,'1H',[-1 2]*1.5 + 0.8);
Spins = nucspinadd(Spins,'1H',[-1 2]*0.5 + 0.7);
Spins.lwEndor = 0.1;

% Experimental parameters
Experiment.Field = 308.46;
Experiment.Range = [6 20];
Experiment.mwFreq = 9.681;

% return each transition separately
Options.Output = 'separate';

[freq,spec_sep] = salt(Spins,Experiment,Options);
spec_sum = sum(spec_sep,1)/4;

% Display
subplot(2,1,1);
plot(freq,spec_sum); axis tight; title('summed');
subplot(2,1,2);
plot(freq,spec_sep); axis tight; title('separate');
xlabel('frequency [MHz]');
