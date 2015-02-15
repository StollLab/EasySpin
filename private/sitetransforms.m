function out = sitetransforms(ID,vec)
% Provides sets of rotation matrices of the rotation point group belonging
% to the space group given in ID. The rotation matrices are used to transform
% tensors between different equivalent sites.

% Input:
%   ID    space group symbol or number, or point group symbol
%   vec   vector to transform (optional)
% Output:
%   out   either cell array of active rotation matrices
%         such that vrot = R{k}*vec is the actively rotated vector vec
%         or an array with vrot for each R{k}, if vec is given

persistent SpaceGroupNames
if isempty(SpaceGroupNames)
  EasySpinPath = fileparts(which(mfilename));
  SpaceGroupDataFile = [EasySpinPath filesep 'spacegroups.txt'];
  [SpaceGroupNo,SpaceGroupNames] = ...
    textread(SpaceGroupDataFile,'%d %s','commentstyle','matlab');
end
PointGroupsSchoenflies = {'C1','Ci','C2','Cs','C2h','D2','C2v','D2h',...
  'C4','S4','C4h','D4','C4v','D2d','D4h','C3','C3i','D3','C3v','D3d',...
  'C6','C3h','C6h','D6','C6v','D3h','D6h','T','Th','O','Td','Oh'};
PointGroupsHermannMauguin = {'1','-1','2','m','2/m','222','mm2','mmm',...
  '4','-4','4/m','422','4mm','-42m','4/mmm','3','-3','32','3m','-3m',...
  '6','-6','6/m','622','6mm','-6m2','6/mmm','23','m-3','432','-43m','m-3m'};
LaueString = {'L1','L2','L3','L4','L5','L6','L7','L8','L9','L10','L11'};
LaueClasses = [1 1 2 2 2 3 3 3 4 4 4 5 5 5 5 6 6 7 7 7 8 8 8 9 9 9 9 10 10 11 11 11];

% Determine Laue class from space group number, space group
% symbol, or point group symbol
%--------------------------------------------------------
LaueClass = [];
if ischar(ID)
  idx = find(strcmp(ID,PointGroupsSchoenflies));
  if isempty(idx)
    idx = find(strcmp(ID,PointGroupsHermannMauguin));
  end
  if ~isempty(idx)
    LaueClass = LaueClasses(idx);
  else
    idx = find(strcmp(ID,SpaceGroupNames));
    if ~isempty(idx)
      ID = idx;
    else
      idx = find(strcmp(ID,LaueString));
      if ~isempty(idx)
        LaueClass = idx;
      else
        error('Unknown point or space group symmetry ''%s''',ID);
      end
    end
  end
end

if isempty(LaueClass)
  if ~isnumeric(ID)
    error('Space group number between 1 and 230.');
  elseif numel(ID)~=1
    error('Space group ID must be a number between 1 and 230.');
  elseif mod(ID,1)~=0 || ~isreal(ID) || (ID<1) || (ID>230)
    error('Wrong space group number %g. Must be between 1 and 230.',ID);
  end
  if     ID<=  2, LaueClass =  1; % (C1, Ci=S2), triclinic
  elseif ID<= 15, LaueClass =  2; % (C2, Cs=C1h, C2h), monoclinic
  elseif ID<= 74, LaueClass =  3; % (D2, C2v, D2h), orthorhombic
  elseif ID<= 88, LaueClass =  4; % (C4, S4, C4h), tetragonal
  elseif ID<=142, LaueClass =  5; % (D4, C4v, D2d, D4h), tetragonal
  elseif ID<=148, LaueClass =  6; % (C3, C3i=S6), trigonal
  elseif ID<=167, LaueClass =  7; % (D3, C3v, D3d), trigonal
  elseif ID<=176, LaueClass =  8; % (C6, C3h, C6h), hexagonal
  elseif ID<=194, LaueClass =  9; % (D6, C6v, D3h, D6h), hexagonal
  elseif ID<=206, LaueClass = 10; % (T, Th), cubic
  elseif ID<=230, LaueClass = 11; % (O, Td, Oh), cubic
  end
end

% Compile list of rotational symmetry transformation for the given Laue class
%-----------------------------------------------------------------------------
% These are matrix representations of the proper rotation group from the given
% Laue class, e.g. D3 for the class (D3, C3v, D3d).
% The proper rotation groups are C1, C2, D2, C4, D4, C3, D3, C6, D6, T, and O.
% (cyclic: C1, C2, C3, C4, C6; dihedral: D2, D3, D4, D6; cubic: T, O)

% For a list of the matrices, see Table II in
%   J.A. Weil, T. Buch, J.E. Clapp, Adv. Magn. Reson. 6, 183-257 (1973)

% All rotations are _active_ rotations: When R is applied to vector v (vrot = R*v),
% then the result vrot is the vector rotated around the specified axis by the
% specified angle.

% Pre-define the most common rotation matrices
E = [+1 0 0; 0 +1 0; 0 0 +1]; % identity operation
C2x = [+1 0 0; 0 -1 0; 0 0 -1]; % C2 around (1,0,0) = x
C2y = [-1 0 0; 0 +1 0; 0 0 -1]; % C2 around (0,1,0) = y
C2z = [-1 0 0; 0 -1 0; 0 0 +1]; % C2 around (0,0,1) = z
C4zp = [0 -1 0; +1 0 0; 0 0 +1]; % C4+ around z
C4zm = [0 +1 0; -1 0 0; 0 0 +1]; % C4- around z
C2d1 = [0 +1 0; +1 0 0; 0 0 -1]; % C2 around (1,+1,0)
C2d2 = [0 -1 0; -1 0 0; 0 0 -1]; % C2 around (1,-1,0)
C3zp = [-1/2 -sqrt(3)/2 0; +sqrt(3)/2 -1/2 0; 0 0 +1]; % C3+ around z
C3zm = [-1/2 +sqrt(3)/2 0; -sqrt(3)/2 -1/2 0; 0 0 +1]; % C3- around z

switch LaueClass
  case 1, % C1 (C1, Ci=S2)
    R{1} = E;
  case 2, % C2 (C2, Cs=C1h, C2h)
    R{1} = E;
    R{2} = C2z;
  case 3, % D2 (D2, C2v, D2h)
    R{1} = E;
    R{2} = C2z;
    R{3} = C2x;
    R{4} = C2y;
  case 4, % C4 (C4, S4, C4h)
    R{1} = E;
    R{2} = C2z;
    R{3} = C4zp;
    R{4} = C4zm;
  case 5, % D4 (D4, C4v, D2d, D4h)
    R{1} = E;
    R{2} = C2z;
    R{3} = C4zp;
    R{4} = C4zm;
    R{5} = C2x;
    R{6} = C2y;
    R{7} = C2d1;
    R{8} = C2d2;
  case 6, % C3 (C3, C3i=S6)
    R{1} = E;
    R{2} = C3zp;
    R{3} = C3zm;
  case 7, % D3 (D3, C3v, D3d)
    R{1} = E;
    R{2} = C3zp;
    R{3} = C3zm;
    R{4} = C2x;
    R{5} = [-1/2 +sqrt(3)/2 0; +sqrt(3)/2 +1/2 0; 0 0 -1];  % C2 (1,+sqrt(3),0)
    R{6} = [-1/2 -sqrt(3)/2 0; -sqrt(3)/2 +1/2 0; 0 0 -1];  % C2 (1,-sqrt(3),0)
  case 8, % C6 (C6, C3h, C6h)
    R{1} = E;
    R{2} = C2z;
    R{3} = C3zp;
    R{4} = C3zm;
    R{5} = [+1/2 -sqrt(3)/2 0; +sqrt(3)/2 +1/2 0; 0 0 +1];  % C6+ z
    R{6} = [+1/2 +sqrt(3)/2 0; -sqrt(3)/2 +1/2 0; 0 0 +1];  % C6- z
  case 9, % D6 (D6, C6v, D3h, D6h)
    R{1} = E;
    R{2} = C2z;
    R{3} = C3zp;
    R{4} = C3zm;
    R{5} = [+1/2 -sqrt(3)/2 0; +sqrt(3)/2 +1/2 0; 0 0 +1];  % C6+ z
    R{6} = [+1/2 +sqrt(3)/2 0; -sqrt(3)/2 +1/2 0; 0 0 +1];  % C6- z
    R{7} = C2x;
    R{8} = C2y;
    R{9} = [-1/2 +sqrt(3)/2 0; +sqrt(3)/2 +1/2 0; 0 0 -1];  % C2 (1,+sqrt(3),0)
    R{10}= [-1/2 -sqrt(3)/2 0; -sqrt(3)/2 +1/2 0; 0 0 -1];  % C2 (1,-sqrt(3),0)
    R{11}= [+1/2 +sqrt(3)/2 0; +sqrt(3)/2 -1/2 0; 0 0 -1];  % C2 (sqrt(3),+1,0)
    R{12}= [+1/2 -sqrt(3)/2 0; -sqrt(3)/2 -1/2 0; 0 0 -1];  % C2 (sqrt(3),-1,0)
  case 10, % T (T, Th)
    R{1} = E;
    R{2} = C2z;
    R{3} = C2x;
    R{4} = C2y;
    R{5} = [0 0 +1; +1 0 0; 0 +1 0];  % C3+ (+1,+1,+1)
    R{6} = [0 +1 0; 0 0 +1; +1 0 0];  % C3- (+1,+1,+1)
    R{7} = [0 0 -1; -1 0 0; 0 +1 0];  % C3+ (+1,-1,-1)
    R{8} = [0 -1 0; 0 0 +1; -1 0 0];  % C3- (+1,-1,-1)
    R{9} = [0 0 +1; -1 0 0; 0 -1 0];  % C3+ (-1,+1,-1)
    R{10}= [0 -1 0; 0 0 -1; +1 0 0];  % C3- (-1,+1,-1)
    R{11}= [0 0 -1; +1 0 0; 0 -1 0];  % C3+ (-1,-1,+1)
    R{12}= [0 +1 0; 0 0 -1; -1 0 0];  % C3- (-1,-1,+1)
  case 11, % O (O, Td, Oh)
    R{1} = E;
    R{2} = C2z;
    R{3} = C4zp;
    R{4} = C4zm;
    R{5} = C2x;
    R{6} = C2y;
    R{7} = C2d1;
    R{8} = C2d2;
    R{9} = [0 0 +1; +1 0 0; 0 +1 0];  % C3+ (+1,+1,+1)
    R{10}= [0 +1 0; 0 0 +1; +1 0 0];  % C3- (+1,+1,+1)
    R{11}= [0 0 -1; -1 0 0; 0 +1 0];  % C3+ (+1,-1,-1)
    R{12}= [0 -1 0; 0 0 +1; -1 0 0];  % C3- (+1,-1,-1)
    R{13}= [0 0 +1; -1 0 0; 0 -1 0];  % C3+ (-1,+1,-1)
    R{14}= [0 -1 0; 0 0 -1; +1 0 0];  % C3- (-1,+1,-1)
    R{15}= [0 0 -1; +1 0 0; 0 -1 0];  % C3+ (-1,-1,+1)
    R{16}= [0 +1 0; 0 0 -1; -1 0 0];  % C3- (-1,-1,+1)
    R{17}= [+1 0 0; 0 0 -1; 0 +1 0];  % C4+ x
    R{18}= [+1 0 0; 0 0 +1; 0 -1 0];  % C4- x
    R{19}= [0 0 +1; 0 +1 0; -1 0 0];  % C4+ y
    R{20}= [0 0 -1; 0 +1 0; +1 0 0];  % C4- y
    R{21}= [0 0 +1; 0 -1 0; +1 0 0];  % C2 (1,0,+1)
    R{22}= [0 0 -1; 0 -1 0; -1 0 0];  % C2 (1,0,-1)
    R{23}= [-1 0 0; 0 0 +1; 0 +1 0];  % C2 (0,1,+1)
    R{24}= [-1 0 0; 0 0 -1; 0 -1 0];  % C2 (0,1,-1)
end

if (nargin==1)
  % Return set of rotation matrices
  out = R;
else
  % Apply site transformations to input vector
  nRotationMatrices = size(R,1);
  for iR = 1:nRotationMatrices
    vecrot(:,iR) = R{iR}*vec(:);
  end
  out = vecrot;
end
