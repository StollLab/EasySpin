% Effect of modulation amplitude on the line shape in cw EPR
%==========================================================================
clear, clf

% Define a single-line absorption spectrum
B = linspace(330,370,1e3);  % field range, mT
lwpp = 2; % line width (peak-to-eak), mT
FWHM = 2/sqrt(3)*lwpp; % line width (full width at half height), mT
spc0 = lorentzian(B,350,FWHM);

% Compute the 1st-harmonic spectrum detected with field modulation
modamp = [0.05 0.1 0.2 0.5 1 2 5]*lwpp;  % peak-to-peak mod.amp., mT
for m = 1:numel(modamp)
  spc = fieldmod(B,spc0,modamp(m),1); % 1st harmonic
  modspc(m,:) = spc;
  modspc_nrm(m,:) = spc/max(spc); % normalized 1st harmonic
end

% Plot results
subplot(3,1,1);
plot(B,spc0);
title('Original absorption signal');
axis tight

subplot(3,1,2);
plot(B,modspc);
title('First-harmonic signal with various modulation amplitudes (unnormalized)');
axis tight

subplot(3,1,3);
plot(B,modspc_nrm);
title('First-harmonic signal with various modulation amplitudes (normalized)');
%axis tight
xlabel('field (mT)')
legend('0.05','0.1','0.2','0.5','1','2','5')

% A Lorentzian line shape is distorted in standard cw EPR
% measurement of the first harmonic if the the peak-to-peak modulation
% amplitude exceeds about 20% of the peak-to-peak linewidth.
