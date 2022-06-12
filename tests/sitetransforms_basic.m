function ok = test()

% Assert that the correct number of rotation matrices is returned
% for every space group.

nSites = zeros(1,230);
for ID = 1:230
  R = runprivate('sitetransforms',ID);
  nSites(ID) = numel(R);
end

nSites_ref = zeros(1,230);
nSites_ref(  1:  2) = 1;  % (C1, Ci=S2), triclinic
nSites_ref(  3: 15) = 2;  % (C2, Cs=C1h, C2h), monoclinic
nSites_ref( 16: 74) = 4;  % (D2, C2v, D2h), orthorhombic
nSites_ref( 75: 88) = 4;  % (C4, S4, C4h), tetragonal
nSites_ref( 89:142) = 8;  % (D4, C4v, D2d, D4h), tetragonal
nSites_ref(143:148) = 3;  % (C3, C3i=S6), trigonal
nSites_ref(149:167) = 6;  % (D3, C3v, D3d), trigonal
nSites_ref(168:176) = 6;  % (C6, C3h, C6h), hexagonal
nSites_ref(177:194) = 12; % (D6, C6v, D3h, D6h), hexagonal
nSites_ref(195:206) = 12; % (T, Th), cubic
nSites_ref(207:230) = 24; % (O, Td, Oh), cubic

ok = nSites==nSites_ref;
