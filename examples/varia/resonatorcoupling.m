% Tuning picture of a cw EPR spectrometer
%===============================================================
% See Krymov, Gerfen, J.Magn.Reson. 162 (2003) 466-478
% https://doi.org/10.1016/S1090-7807(03)00109-5
% mainly Eq. (7)

clear, clf

% Coupling coefficient beta
% (1 critical coupling, <1 undercoupled, >1 overcoupled)
beta = [0.2 1 3 8];

% Scaled frequency offset xi
xi = linspace(-1,1,301)*20;

% Compute complex-valued reflection coefficient Gamma
for b = 1:numel(beta)
  Gamma(b,:) = (1-beta(b)-1i*xi)./(1+beta(b)-1i*xi); 
end

% Plotting
plot(xi,abs(Gamma));
xlabel('scaled frequency offset');
ylabel('reflection');
legend('undercoupled','critically coupled','overcoupled','overcoupled');
