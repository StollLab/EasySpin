% mhz2mt   obsolete
% use unitconvert(input, 'MHz->mT', g)

function varargout = mhz2mt(x_MHz,g)

if nargin<1
  error('The function mhz2mt is obsolete. Use:\noutput = unitconvert(1, ''MHz->mT'')');
end
if nargin<2
  error('The function mhz2mt is obsolete. Use:\noutput = unitconvert(%d, ''MHz->mT'')',x_MHz);
end

error('The function mhz2mt is obsolete. Use:\noutput = unitconvert(%d, ''MHz->mT'', %d)',x_MHz,g);

end
