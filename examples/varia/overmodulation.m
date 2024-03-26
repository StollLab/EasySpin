% Effect of modulation amplitude on the line shape in cw EPR
%==========================================================================
clear, clf

% Define a single-line absorption spectrum
B = linspace(300,400,1e3);  % mT
FWHM = 5; % line width (full width at half height)
spc0 = lorentzian(B,350,FWHM);

% Compute the 1st-harmonic signal detected with field modulation
modamp = [0.05 0.1 0.2 0.5 1 2 5]*FWHM;
for m = 1:numel(modamp)
  spc = fieldmod(B,spc0,modamp(m),1); % 1st harmonic
  modspc(m,:) = spc/max(spc); % normalized 1st harmonic
end

% Plot results
subplot(2,1,1);
plot(B,spc0);
title('Original absorption signal');
axis tight

subplot(2,1,2);
plot(B,modspc);
title('1st harmonic signal with various modulation amplitudes');
axis tight

% Conclusion: A Lorentzian lnie shape is distorted in
% standard cw EPR measurement of the first harmonic,
% if the the peak-to-peak modulation amplitude exceeds
% 20% of the FWHM of the absorption line.
