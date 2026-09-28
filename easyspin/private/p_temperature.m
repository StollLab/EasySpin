% p_temperature  Validate Exp.Temperature
%
%   useTemperature = p_temperature(Exp)
%
%   Input:
%     Exp             experiment structure
%
%   Output:
%     useTemperature  false if Exp.Temperature is absent, empty, or NaN
%                     (high-temperature limit), true otherwise
%
%   If given, Exp.Temperature must be a single real, finite, non-negative
%   number (in K). Otherwise, an error is issued.

function useTemperature = p_temperature(Exp)

useTemperature = false;

if ~isfield(Exp,'Temperature') || isempty(Exp.Temperature)
  return
end

T = Exp.Temperature;
if isnumeric(T) && isscalar(T) && isnan(T)
  return
end

if ~isnumeric(T) || ~isscalar(T) || ~isreal(T) || ~isfinite(T) || T<0
  error('Exp.Temperature must be a single non-negative finite number (in K). Omit it or set it to NaN for the high-temperature limit.');
end

useTemperature = true;

end
