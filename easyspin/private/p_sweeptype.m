% p_sweeptype  Determine sweep type and validate Exp.mwFreq, Exp.Field, Exp.nPoints
%
%   FieldSweep = p_sweeptype(Exp,allowZeroField)
%
%   Input:
%     Exp             experiment structure
%     allowZeroField  true if frequency sweeps at Exp.Field=0 are supported
%
%   Output:
%     FieldSweep      true for field sweeps (Exp.mwFreq given), false for
%                     frequency sweeps (Exp.Field given)
%
%   Also rejects obsolete fields and sweep range fields that don't match
%   the sweep type (e.g. Exp.mwRange in a field sweep). Empty or NaN range
%   fields count as absent.

function FieldSweep = p_sweeptype(Exp,allowZeroField)

% Obsolete fields
if isfield(Exp,'Orientations')
  error('Exp.Orientations is no longer supported, use Exp.SampleFrame instead.');
end
if isfield(Exp,'CrystalOrientation')
  error('Exp.CrystalOrientation is no longer supported, use Exp.SampleFrame/Exp.MolFrame instead.');
end
if isfield(Exp,'Mode')
  error('Exp.Mode is no longer supported. Use Exp.mwMode instead.');
end

% Field sweep (Exp.mwFreq) or frequency sweep (Exp.Field)
mwFreqGiven = isfield(Exp,'mwFreq') && ~isempty(Exp.mwFreq);
FieldGiven = isfield(Exp,'Field') && ~isempty(Exp.Field);
if mwFreqGiven==FieldGiven
  msg = 'Provide either Exp.mwFreq (for a field sweep) or Exp.Field (for a frequency sweep)';
  if mwFreqGiven
    error('%s, but not both.',msg);
  else
    error('%s.',msg);
  end
end
FieldSweep = mwFreqGiven;

% Sweep range fields must match the sweep type
if FieldSweep
  wrongFields = {'mwRange','mwCenterSweep'};
  sweepStr = 'a field sweep (Exp.mwFreq is given)';
else
  wrongFields = {'Range','CenterSweep'};
  sweepStr = 'a frequency sweep (Exp.Field is given)';
end
for k = 1:numel(wrongFields)
  f = wrongFields{k};
  if isfield(Exp,f) && ~isempty(Exp.(f)) && ~all(isnan(Exp.(f)))
    error('Exp.%s cannot be used in %s.',f,sweepStr);
  end
end

% Microwave frequency or static field
if FieldSweep
  mwFreq = Exp.mwFreq;
  if ~isnumeric(mwFreq) || ~isscalar(mwFreq) || ~isreal(mwFreq) || ~isfinite(mwFreq) || mwFreq<=0
    error('Uninterpretable microwave frequency in Exp.mwFreq.');
  end
  logmsg(1,'  field sweep, mw frequency %0.8g GHz',mwFreq);
else
  Field = Exp.Field;
  if ~isnumeric(Field) || ~isscalar(Field) || ~isreal(Field) || ~isfinite(Field)
    error('Uninterpretable magnetic field in Exp.Field.');
  end
  if Field<0
    error('Exp.Field cannot be negative. Negative fields are only supported for field sweeps.');
  end
  if Field==0 && ~allowZeroField
    error('Exp.Field cannot be zero.');
  end
  logmsg(1,'  frequency sweep, magnetic field %0.8g mT',Field);
end

% Number of points
nPoints = Exp.nPoints;
if ~isnumeric(nPoints) || ~isscalar(nPoints) || ~isreal(nPoints) || ~isfinite(nPoints) || ...
    nPoints~=round(nPoints) || nPoints<2
  error('Problem with Exp.nPoints. Needs to be an integer not smaller than 2.');
end

end
