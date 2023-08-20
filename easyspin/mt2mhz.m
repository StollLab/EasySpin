% mt2mhz   obsolete
% use unitconvert(input, 'mT->MHz')

function mt2mhz(x_mT,g)

if nargin<1, error(sprintf('The function mhz2mt is obsolete. Use:\noutput = unitconvert(1, ''mT->MHz'')')); end
if nargin<2
    error(sprintf('The function mt2mhz is obsolete. Use:\noutput = unitconvert(%d, ''mT->MHz'')',x_mT));
end

error(sprintf('The function mt2mhz is obsolete. Use:\noutput = unitconvert(%d, ''mT->MHz'', %d)',x_mT,g));

end