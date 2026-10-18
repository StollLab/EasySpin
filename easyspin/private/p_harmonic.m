% p_harmonic  Validate Exp.Harmonic and Exp.ModAmp
%
%   [Exp,derivsLeft] = p_harmonic(Exp,FieldSweep,hasBroadening)
%
%   Input:
%     Exp            experiment structure
%     FieldSweep     true for field sweeps, false for frequency sweeps
%     hasBroadening  true if the spectrum has any line broadening (field
%                    modulation via Exp.ModAmp counts as broadening, too)
%
%   Output:
%     Exp            experiment structure, with Harmonic and ModAmp set
%     derivsLeft     number of derivatives still to be applied to the
%                    absorption spectrum (0 with field modulation, since
%                    fieldmod generates the harmonic)
%
%   If Exp.Harmonic is absent, empty or NaN, it defaults to 1 for field
%   sweeps with broadening and to 0 otherwise.

function [Exp,derivsLeft] = p_harmonic(Exp,FieldSweep,hasBroadening)

% Field modulation amplitude
if ~isfield(Exp,'ModAmp') || isempty(Exp.ModAmp)
  Exp.ModAmp = 0;
end
if ~isnumeric(Exp.ModAmp) || ~isscalar(Exp.ModAmp) || isnan(Exp.ModAmp) || Exp.ModAmp<0
  error('Exp.ModAmp must be either a single positive number or zero.');
end
fieldModulation = Exp.ModAmp>0;
if fieldModulation && ~FieldSweep
  error('Exp.ModAmp cannot be used with frequency sweeps.');
end

% Detection harmonic
if ~isfield(Exp,'Harmonic') || isempty(Exp.Harmonic) || any(isnan(Exp.Harmonic))
  if FieldSweep && (hasBroadening || fieldModulation)
    Exp.Harmonic = 1;
  else
    Exp.Harmonic = 0;
  end
end
if ~isnumeric(Exp.Harmonic) || ~isscalar(Exp.Harmonic) || ~any(Exp.Harmonic==[0 1 2])
  error('Exp.Harmonic must be 0, 1 or 2.');
end

if fieldModulation
  if Exp.Harmonic<1
    error('With field modulation (Exp.ModAmp), Exp.Harmonic=0 does not work.');
  end
  derivsLeft = 0;
  logmsg(1,'  harmonic %d, field modulation with amplitude %g mT',Exp.Harmonic,Exp.ModAmp);
else
  if ~hasBroadening && Exp.Harmonic~=0
    error('No broadening given. Cannot compute spectrum with Exp.Harmonic=%d.',Exp.Harmonic);
  end
  derivsLeft = Exp.Harmonic;
  logmsg(1,'  harmonic %d',Exp.Harmonic);
end

end
