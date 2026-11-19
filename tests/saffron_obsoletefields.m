function ok = test()  %#ok<FNDEF>

% The old fields for the pulse sequence (Exp.Flip, Exp.Phase, Exp.tp, Exp.t,
% Exp.Inc) are obsolete: saffron rejects them for predefined and
% user-defined experiments

Sys.S = 1/2;
Sys.g = 2;

ExpPredefined.Field = 350;  % mT
ExpPredefined.mwFreq = 9.8;  % GHz
ExpPredefined.Sequence = '2pESEEM';
ExpPredefined.tau = 0.1;  % µs
ExpPredefined.dt = 0.02;  % µs
ExpPredefined.nPoints = 3;

p90.Flip = pi/2;
p180.Flip = pi;
ExpCustom = rmfield(ExpPredefined,{'tau','dt'});
ExpCustom.Sequence = {p90 0.1 p180 0.1};
ExpCustom.Dim1 = {'d1,d2' 0.02};

obsoleteFields = {'Flip','Phase','tp','t','Inc'};
for iField = 1:numel(obsoleteFields)
  field_ = obsoleteFields{iField};
  ExpP = ExpPredefined;
  ExpP.(field_) = [0 0];
  ExpC = ExpCustom;
  ExpC.(field_) = [0 0];
  calls = {@() saffron(Sys,ExpP), @() saffron(Sys,ExpC)};
  for k = 1:numel(calls)
    if strcmp(field_,'tp') && k==1
      % Exp.tp is a parameter of saffron's predefined experiments
      ok(iField,k) = true;
      continue
    end
    try
      calls{k}();
      ok(iField,k) = false;
    catch err
      ok(iField,k) = contains(err.message,['Exp.' field_ ' is no longer supported']);
    end
  end
end
ok = ok(:).';
