% unitconvert   Unit conversion tool
%
%   output = unitconvert(input, units)
%   output = unitconvert(input, units, g)
%
%   Converts the value in input, through the conversion specified in 
%   units returning the result in output.
%
%   units take the form of: 'unit_a->unit_b' where unit_a and unit_b are:
%   cm^-1, eV, K, mT or MHz
%
%   e.g. 'cm^-1->MHz' for conversion of wavenumbers to megahertz 
%
%                'cm^-1->eV'  'cm^-1->K'  'cm^-1->mT'  'cm^-1->MHz'       
%   'eV->cm^-1'               'eV->K'     'eV->mT'     'eV->MHz'
%   'K->cm^-1'   'K->eV'                  'K->mT'      'K->MHz'
%   'mT->cm^-1'  'mT->eV'     'mT->K'                  'mT->MHz'
%   'MHz->cm^-1' 'MHz->eV'    'MHz->K'    'MHz->mT'
%
%   When converting into or from magnetic field units, a g factor given 
%   as the third parameter is used. If it is not given, the g factor
%   of the free electron (gfree) is used.
%
%   input can be a vector of values. In this case, g
%   can be a scalar or a vector of the same size as input.
%
%   Example:
%     value_MHz = unitconvert(value_wn,'cm^-1->MHz')
%     value_mT = unitconvert(value_wn,'cm^-1->mT',2.005)

function out = unitconvert(value,units,g)
if nargin==0 && nargout==0, help(mfilename); return; end
if nargin<=1 && nargout==1
    error('Wrong number of input arguments');
end

if nargin<3, g = gfree; end

% Cell array with unit conversion as string and function handles to compute
% the conversion
% From Matlab 2022b, dictionaries would be better suited

unitconv = {
  "cm^-1->eV",    @(value) value.*100*clight*planck/evolt;
  "cm^-1->K",     @(value) value.*100*clight*planck/boltzm;
  "cm^-1->mT",    @(value) value./g*(planck/bmagn/1e-3)*100*clight;
  "cm^-1->MHz",   @(value) value.*100*clight/1e6;

  "eV->cm^-1",    @(value) value.*evolt/100/clight/planck;
  "eV->K",        @(value) value.*evolt/boltzm;
  "eV->mT",       @(value) value./g/bmagn/1e-3*evolt;
  "eV->MHz",      @(value) value.*evolt/planck/1e6;

  "K->cm^-1",     @(value) value.*boltzm/100/clight/planck;
  "K->eV",        @(value) value.*boltzm/evolt;
  "K->mT",        @(value) value./g/bmagn/1e-3*boltzm;
  "K->MHz",       @(value) value.*boltzm/planck/1e6;

  "mT->cm^-1",    @(value) value.*g/(planck/bmagn/1e-3)/100/clight;
  "mT->eV",       @(value) value.*g*bmagn*1e-3/evolt;
  "mT->K",        @(value) value.*g*bmagn*1e-3/boltzm;
  "mT->MHz",      @(value) value.*g*(1e-3*bmagn/planck/1e6);

  "MHz->cm^-1",   @(value) value.*1e6/100/clight;
  "MHz->eV",      @(value) value.*1e6*planck/evolt;
  "MHz->K",       @(value) value.*1e6*planck/boltzm;
  "MHz->mT",      @(value) value./g*(planck/bmagn/1e-3)*1e6
  };

% List of possible unit conversion as string array
unames = string(unitconv(:,1));

% Test if input unit conversion exists (case-sensitive)
unitMatch = strcmp(unames,units);

if any(unitMatch)
  % Return unit conversion
  out = unitconv{unitMatch,2}(value);
else
  % Throw errors depending on unit input
  % Test if input unit conversion exists (case-insensitive)
  unitMatchI = strcmpi(unames,units);
  if any(unitMatchI)
    error('You provided: %s. Did you mean %s?',units,unitconv{unitMatchI,1})
  else
    error('Unknown unit conversion specified.')
  end
end