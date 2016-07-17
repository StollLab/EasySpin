% ENDOR orientation selection in a powder
%=====================================================================
clear, clf

% Spin system
Sys.S = 1/2;
Sys.g = [2.3 2.1 2];
Sys = nucspinadd(Sys,'1H',[3,6,2]);
Sys.HStrain = [1 1 1]*220;     % MHz
Sys.lwEndor = 0.2;    % MHz

% ENDOR experiment settings
Exp.mwFreq = 9.5;
Exp.ExciteWidth = 100;
Exp.Range = [8 20];

Fields = 295:5:340;

for iField = 1:numel(Fields)   % loop over all field values
  Exp.Field = Fields(iField);
  [freq,spectra(iField,:)] = salt(Sys,Exp);
end

% Stack plot
stackplot(freq,spectra);
xlabel('frequency (MHz)');
title('ENDOR spectra');

