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

switch lower(units)

    case 'cm^-1->ev'
        out = value.*100*clight*planck/evolt;
    case 'cm^-1->k'
        out = value.*100*clight*planck/boltzm;
    case 'cm^-1->mt'
        out = value./g*(planck/bmagn/1e-3)*100*clight;
    case 'cm^-1->mhz'
        out = value.*100*clight/1e6;


    case 'ev->cm^-1'
        out = value.*evolt/100/clight/planck;
    case 'ev->k'
        out = value.*evolt/boltzm;
    case 'ev->mt'
        out = value./g/bmagn/1e-3*evolt;
    case 'ev->mhz'
        out = value.*evolt/planck/1e6;


    case 'k->cm^-1'
        out = value.*boltzm/100/clight/planck;
    case 'k->ev'
        out = value.*boltzm/evolt;
    case 'k->mt'
        out = value./g/bmagn/1e-3*boltzm;
    case 'k->mhz'
        out = value.*boltzm/planck/1e6;


    case 'mt->cm^-1'
        out = value.*g/(planck/bmagn/1e-3)/100/clight;
    case 'mt->ev'
        out = value.*g*bmagn*1e-3/evolt;
    case 'mt->k'
        out = value.*g*bmagn*1e-3/boltzm;
    case 'mt->mhz'
        out = value.*g*(1e-3*bmagn/planck/1e6);


    case 'mhz->cm^-1'
        out = value.*1e6/100/clight;
    case 'mhz->ev'
        out = value.*1e6*planck/evolt;
    case 'mhz->k'
        out = value.*1e6*planck/boltzm;
    case 'mhz->mt'
        out = value./g*(planck/bmagn/1e-3)*1e6;
        

    otherwise
        error('Unknown unit conversion specified ')

end
