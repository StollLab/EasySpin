% p_sweeprange  Sweep range from Exp.CenterSweep/Exp.Range or Exp.mwCenterSweep/Exp.mwRange
%
%   Range = p_sweeprange(Exp,freqsweep,allowNegative)
%
%   Input:
%     Exp            experiment structure
%     freqsweep      false: use Exp.CenterSweep and Exp.Range
%                    true:  use Exp.mwCenterSweep and Exp.mwRange
%     allowNegative  true if negative values are allowed in the range
%
%   Output:
%     Range          sweep range [lo hi], or [] if neither field is given
%
%   Exp.CenterSweep (Exp.mwCenterSweep) has precedence over Exp.Range
%   (Exp.mwRange). A field is considered not given if it is absent, empty,
%   or NaN. An error is issued if the range is invalid, or if it contains
%   negative values and allowNegative is false.

function Range = p_sweeprange(Exp,freqsweep,allowNegative)

if freqsweep
  csField = 'mwCenterSweep';
  rangeField = 'mwRange';
else
  csField = 'CenterSweep';
  rangeField = 'Range';
end

isGiven = @(f) isfield(Exp,f) && ~isempty(Exp.(f)) && ~all(isnan(Exp.(f)(:)));

if isGiven(csField)
  CenterSweep = Exp.(csField);
  if ~isnumeric(CenterSweep) || numel(CenterSweep)~=2
    error('Invalid sweep range! Check Exp.%s or Exp.%s.',csField,rangeField);
  end
  Range = CenterSweep(1) + [-1 1]*CenterSweep(2)/2;
elseif isGiven(rangeField)
  Range = Exp.(rangeField);
  if ~isnumeric(Range) || numel(Range)~=2
    error('Invalid sweep range! Check Exp.%s or Exp.%s.',csField,rangeField);
  end
  Range = reshape(Range,1,2);
else
  Range = [];
  return
end

if ~isreal(Range) || any(~isfinite(Range)) || Range(1)>=Range(2)
  error('Invalid sweep range! Check Exp.%s or Exp.%s.',csField,rangeField);
end

if ~allowNegative && any(Range<0)
  error('Sweep range cannot be negative! Check Exp.%s or Exp.%s.',csField,rangeField);
end

end
