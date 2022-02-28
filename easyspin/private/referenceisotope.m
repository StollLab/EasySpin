% Determine hyperfine and quadrupole reference isotopes for a given
% element.
%
% The reference isotope for hyperfine coupling A is the most (naturally)
% abundant isotope with spin 1/2 or larger.
% The reference isotope for quadrupole coupling Q is the most (naturally)
% abundant isotope with spin 1 or larger.

function [gref,qref] = referenceisotope(el)

% Convert element symbol (such as "Cu" or 'Cu') to element number (29)
if ischar(el) || isstring(el)
  no = elementsymbol2no(el);
else
  no = el;
end

% Get data for all isotopes of given element
data = nucdata;
idx = find(data.Protons==no);
abundances = data.Abundances(idx);
I = data.Spins(idx);
gn = data.gns(idx);
qm = data.qms(idx);
symb = data.Symbols(idx);

% gn reference isotope: most (naturally) abundant isotope with I>=1/2
ab = abundances;
ab(I<1/2) = 0;
[mx,idxg] = max(ab);
if mx>0
  gref.gn = gn(idxg);
  gref.symbol = symb{idxg};
  gref.I = I(idxg);
else
  gref = [];
end

% Q reference isotope: most (naturally) abundant isotope with I>=1
ab = abundances;
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
