% Tuning picture of a cw EPR spectrometer
%===============================================================
% See Krymov, Gerfen, J.Magn.Reson. 162 (2003) 466-478
% mainly Eq. (7)

clear, clf

% coupling coefficient
% (1 critical coupling, <1 undercoupled, >1 overcoupled)
beta = [0.2 1 3 8];

% scaled frequency offset
xi = linspace(-1,1,301)*20;

% compute complex reflection coefficient
for b = 1:numel(beta)
  Gamma(b,:) = (1-beta(b)-1i*xi)./(1+beta(b)-1i*xi); 
end

% display
plot(xi,abs(Gamma));
xlabel('scaled frequency offset');
ylabel('reflection');
legend('undercoupled','critically coupled','overcoupled','overcoupled');
%legend boxoff
