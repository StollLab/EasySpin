% Various methods for de-noising spectral data
%==========================================================================
% Flat and binomial moving averages, Savitzky-Golay and RC filtering.

clear, clf

% Construct a noisy signal.
%--------------------------------------------------------------------------
x = linspace(-1,1,801);
y_pure = gaussian(x,0,0.2);
y_pure = y_pure/max(y_pure);
y = y_pure + 0.1*(rand(size(x))-.5);

% Apply various smoothing methods.
%--------------------------------------------------------------------------
m = 50; % half-width of smoothing window
y_flat = datasmooth(y,m,'flat'); % flat moving average
y_binom = datasmooth(y,m,'binom'); % binomial moving average
y_sg = datasmooth(y,m,'savgol',2); % Savitzky-Golay moving average
y_rc = rcfilt(y,1,m); % RC filter

% Display.
%--------------------------------------------------------------------------
% In the resulting plot it can be seen, that RC filtering
% and standard moving average distort width and amplitude
% of the signal, whereas binomial and Savitzky-Golay
% averages work best.

plot(x,[y+0.2; y_flat; y_binom; y_sg; y_rc]);
legend('none','flat','binomial','Savitzky-Golay','RC');
title('Various smoothing methods');
