% Determine hyperfine and quadrupole reference isotopes for a given
% element.
%
function [gref,qref] = referenceisotope(el)

if ischar(el)
  no = elementsymbol2no(el);
else
  no = el;
end

data = nucdata;
idx = find(data.Protons==no);

ab = data.Abundances(idx);
I = data.Spins(idx);
gn = data.gns(idx);
qm = data.qms(idx);
symb = data.Symbols(idx);

% gn reference isotope: most abundant I>=1/2 isotope
ab(I<1/2) = 0;
[mx,idxg] = max(ab);
if mx>0
  gref.gn = gn(idxg);
  gref.symbol = symb{idxg};
  gref.I = I(idxg);
else
  gref = [];
end

% Q reference isotope: most abundant I>=1 isotope
ab(I<1) = 0;
[mx,idxq] = max(ab);
if mx>0
  qref.qm = qm(idxq);
  qref.symbol = symb{idxq};
  qref.I = I(idxq);
else
  qref = [];
end

end
