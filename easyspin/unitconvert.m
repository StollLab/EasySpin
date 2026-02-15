% unitconvert   Energy and field unit conversion
%
%   vals_out = unitconvert(vals, fromto)
%   vals_out = unitconvert(vals, fromto, g)
%
%   Converts values between energy and field units.
%
%   Input:
%     vals      values to convert, scalar or array
%     fromto    source and target unit, in the form 'unit1->units2', for
%               example 'cm^-1->GHz'. Supported units are:
%                  cm^-1, eV, K
%                  MHz, GHz, THz
%                  G, mT, T
%     g         g factor for conversion from and to magnetic-field units;
%               set to gfree if not provided. If vals is an array, g can
%               be a scalar or an array of the same size as vals
%
%  Output:
%     vals_out  values converted from source to target unit
%
%  Example:
%    value_GHz = unitconvert(value_wn,'cm^-1->GHz')
%    value_mT = unitconvert(value_wn,'cm^-1->mT',2.005)

function value_out = unitconvert(value,units,g)

if nargin==0
  help(mfilename);
  return;
end

if nargin<2
  error('At least two input arguments (values and units) are required.');
end

if nargin<3
  g = gfree;
else
  if ~isnumeric(g)
    error('g value (third input) must be numeric.')
  end
  if numel(g)>1 && numel(value)>1 && ~isequal(size(value),size(g))
    error('If g (third input) is an array, it must be the same size as values (first input).');
  end
end

if ~contains(units,'->')
   error('Unit conversion string must be of the form "unit1->units". "->" is missing.');
end

sourceUnit = extractBefore(units,'->');
targetUnit = extractAfter(units,'->');

supportedUnits = {'cm^-1', 'eV', 'K', 'MHz', 'GHz', 'THz', 'G', 'mT', 'T'};

% Convert values from source unit to canonical intermediate unit (J)
switch sourceUnit
  case 'cm^-1'
    value_J = planck*clight*value*100; 
  case 'eV'
    value_J = evolt*value;
  case 'J'
    value_J = value;
  case 'K'
    value_J = boltzm*value;
  case 'MHz'
    value_J = planck*value*1e6;
  case 'GHz'
    value_J = planck*value*1e9;
  case 'THz'
    value_J = planck*value*1e12;
  case 'G'
    value_J = bmagn*g.*value*1e-4;
  case 'mT'
    value_J = bmagn*g.*value*1e-3;
  case 'T'
    value_J = bmagn*g.*value;
  otherwise
    idx = find(strcmpi(sourceUnit,supportedUnits),1);
    if ~isempty(idx)
      error('Supplied source unit is "%s". Did you mean "%s"?',sourceUnit,supportedUnits{idx});
    else
      error('Unknown source unit "%s".',sourceUnit);
    end
end

% Convert values from canonical intermediate unit (J) to target unit
switch targetUnit
  case 'cm^-1'
    value_out = value_J/planck/clight/100;
  case 'eV'
    value_out = value_J/evolt;
  case 'J'
    value_out = value_J;
  case 'K'
    value_out = value_J/boltzm;
  case 'G'
    value_out = value_J/bmagn./g/1e-4;
  case 'mT'
    value_out = value_J/bmagn./g/1e-3;
  case 'T'
    value_out = value_J/bmagn./g;
  case 'MHz'
    value_out = value_J/planck/1e6;
  case 'GHz'
    value_out = value_J/planck/1e9;
  case 'THz'
    value_out = value_J/planck/1e12;
  otherwise
    idx = find(strcmpi(targetUnit,supportedUnits),1);
    if ~isempty(idx)
      error('Supplied target unit is "%s". Did you mean "%s"?',targetUnit,supportedUnits{idx});
    else
      error('Unknown target unit "%s".',targetUnit);
    end
end

end
