% Effect of modulation amplitude on the line shape in cw EPR
%==========================================================================
clear, clf

% Define a single-line absorption spectrum
x = linspace(300,400,1e3);
FWHM = 5; % line width (full width at half height)
Spc = lorentzian(x,350,FWHM);

% Compute the 1st-harmonic signal detected with field modulation
ModAmpl = [0.05 0.1 0.2 0.5 1 2 5]*FWHM;
for m = 1:numel(ModAmpl)
  y = fieldmod(x,Spc,ModAmpl(m),1); % 1st harmonic
  ModSpc(m,:) = y/max(y); % normalized 1st harmonic
end

% Plot results
subplot(2,1,1);
plot(x,Spc);
title('Original absorption signal');
axis tight

subplot(2,1,2);
plot(x,ModSpc);
title('1st harmonic signal with various modulation amplitudes');
axis tight

% Conclusion: A Lorentzian lnie shape is distorted in
% standard cw EPR measurement of the first harmonic,
% if the the peak-to-peak modulation amplitude exceeds
% 20% of the FWHM of the absorption line.
