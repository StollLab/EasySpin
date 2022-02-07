% sitetransforms   Rotation matrices for space groups
%
%   R = sitetransforms(ID)
%   vrot = sitetransforms(ID,vec)
%
% Provides sets of rotation matrices of the rotation point group belonging
% to the space group given in ID. The rotation matrices are used to
% transform vectors and tensors between different equivalent sites.
%
% Input:
%   ID    One of the following:
%         - Schoenflies or Hermann-Mauguin symbol for crystallographic point group
%         - Hermann-Mauguin symbol for space group (e.g. P121)
%         - Space group number (between 1 and 230)
%   vec   vector to transform (optional)
% Output:
%   R     cell array of active rotation matrices
%         such that vrot = R{k}*vec is the actively rotated vector vec
%   vrot  array with vrot for each R{k}

function out = sitetransforms(ID,vec)

if nargin==0, help(mfilename); return; end

transformVector = nargin==2;
if transformVector
  if numel(vec)~=3
    error('Vector (2nd argument) must have 3 elements.');
  end
end

% Load space group data
%--------------------------------------------------------------------------
persistent SpaceGroups
if isempty(SpaceGroups)
  EasySpinPath = fileparts(which(mfilename));
  SpaceGroupDataFile = [EasySpinPath filesep 'spacegroups.txt'];
  
  % Read in file with space group numbers, symbols and notes.
  % Notes contain information on unique axis and settings.
  fh = fopen(SpaceGroupDataFile);
  C = textscan(fh,'%f %s %s','commentstyle','%');
  fclose(fh);
  [SpaceGroups.No,SpaceGroups.Names,SpaceGroups.Axes] = C{:};
end

% Determine Laue class from input ID
%--------------------------------------------------------------------------

% Process input if it is a character array or string
if ischar(ID)
  idx = find(strcmp(ID,SpaceGroups.Names));
  if ~isempty(idx) % if it's a space group symbol
    LaueClass = groupnumber2laueclass(SpaceGroups.No(idx));
    AxisConvention = SpaceGroups.Axes{idx};
  else
    PointGroupsSchoenflies = {'C1','Ci','C2','Cs','C2h','D2','C2v','D2h',...
      'C4','S4','C4h','D4','C4v','D2d','D4h','C3','C3i','D3','C3v','D3d',...
      'C6','C3h','C6h','D6','C6v','D3h','D6h','T','Th','O','Td','Oh'};
    PointGroupsHermannMauguin = {'1','-1','2','m','2/m','222','mm2','mmm',...
      '4','-4','4/m','422','4mm','-42m','4/mmm','3','-3','32','3m','-3m',...
      '6','-6','6/m','622','6mm','-6m2','6/mmm','23','m-3','432','-43m','m-3m'};
    LaueClasses = [1 1 2 2 2 3 3 3 4 4 4 5 5 5 5 6 6 7 7 7 8 8 8 9 9 9 9 10 10 11 11 11];
    idx = find(strcmp(ID,PointGroupsSchoenflies));
    if isempty(idx)
      idx = find(strcmp(ID,PointGroupsHermannMauguin));
    end
    if ~isempty(idx) % if it's a point group symbol (Schoenflies or H-M)
      LaueClass = LaueClasses(idx);
      AxisConvention = 'default';
    else
      error('Point or space group symmetry symbol ''%s'' is unknown.',ID);
    end
  end
end

% Process input if it is a number
if ~ischar(ID)
  if ~isnumeric(ID)
    error('Space group ID must be a number or a string/character array.');
  elseif numel(ID)~=1
    error('Space group ID must be a single number, not an array.');
  elseif mod(ID,1)~=0 || ~isreal(ID) || ID<1 || ID>230
    error('Invalid space group number %g. Must be between 1 and 230.',ID);
  end
  LaueClass = groupnumber2laueclass(ID);
  AxisConvention = 'default';
end

% Supplement default unique axis
if strcmp(AxisConvention,'default')
  if LaueClass==2 % monoclinic space groups
    AxisConvention = 'b'; % b is commonly used as default -> yC is unique axis
  else
    AxisConvention = 'z'; % set zC as unique axis
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

% Below we also list explicitly how EasySpin's crystal frame axes system
% xC, yC, zC is defined relative to the symmetry axes of the various point
% groups.

% Pre-define the most common rotation matrices
E = [+1 0 0; 0 +1 0; 0 0 +1]; % identity operation
C2x = [+1 0 0; 0 -1 0; 0 0 -1]; % C2 around (1,0,0) = xC
C2y = [-1 0 0; 0 +1 0; 0 0 -1]; % C2 around (0,1,0) = yC
C2z = [-1 0 0; 0 -1 0; 0 0 +1]; % C2 around (0,0,1) = zC
C4zp = [0 -1 0; +1 0 0; 0 0 +1]; % C4+ around zC
C4zm = [0 +1 0; -1 0 0; 0 0 +1]; % C4- around zC
C2xy1 = [0 +1 0; +1 0 0; 0 0 -1]; % C2 around (1,+1,0)
C2xy2 = [0 -1 0; -1 0 0; 0 0 -1]; % C2 around (1,-1,0)
C3zp = [-1/2 -sqrt(3)/2 0; +sqrt(3)/2 -1/2 0; 0 0 +1]; % C3+ around zC
C3zm = [-1/2 +sqrt(3)/2 0; -sqrt(3)/2 -1/2 0; 0 0 +1]; % C3- around zC

switch LaueClass
  case 1 % #1-2, triclinic, C1 (C1, Ci=S2)
    % Axis convention: None. xC, yC, zC are arbitrary.
    R{1} = E;
  case 2 % #3-15, monoclinic, C2 (C2, Cs=C1h, C2h)
    % Axis conventions:
    % (1) Point group given: zC along unique two-fold axis
    % (2) Short HM space group symbol given: yC along two-fold axis
    % (3) Full HM space group symbol given: two-fold axis according to symbol
    R{1} = E;
    switch AxisConvention
      case {'b','-b'},     R{2} = C2y;
      case {'a','-a'},     R{2} = C2x;
      case {'c','-c','z'}, R{2} = C2z;
      otherwise
        error('Unknown unique axis for this monoclinic space group or crystallographic point group.');
    end
  case 3 % #16-74, orthorhombic, D2 (D2, C2v, D2h)
    % Axis conventions: xC, yC, zC along two-fold axes
    R{1} = E;
    R{2} = C2z;
    R{3} = C2x;
    R{4} = C2y;
  case 4 % #75-88, tetragonal, C4 (C4, S4, C4h)
    % Axis convention: zC along four-fold axis, xC and yC arbitrary
    R{1} = E;
    R{2} = C2z;
    R{3} = C4zp;
    R{4} = C4zm;
  case 5 % #89-142, tetragonal, D4 (D4, C4v, D2d, D4h)
    % Axis convention: zC along four-fold axis, xC along one of the two-fold axes
    R{1} = E;
    R{2} = C2z;
    R{3} = C4zp;
    R{4} = C4zm;
    R{5} = C2x;
    R{6} = C2y;
    R{7} = C2xy1;
    R{8} = C2xy2;
  case 6 % #143-148, trigonal, C3 (C3, C3i=S6)
    % Axis convention: zC along three-fold axis, xC and yC arbitrary
    R{1} = E;
    R{2} = C3zp;
    R{3} = C3zm;
  case 7 % #149-167, trigonal, D3 (D3, C3v, D3d)
    % Axis convention: zC along three-fold axis, xC along one of the two-fold axes
    R{1} = E;
    R{2} = C3zp;
    R{3} = C3zm;
    R{4} = C2x;
    R{5} = [-1/2 +sqrt(3)/2 0; +sqrt(3)/2 +1/2 0; 0 0 -1];  % C2 (1,+sqrt(3),0)
    R{6} = [-1/2 -sqrt(3)/2 0; -sqrt(3)/2 +1/2 0; 0 0 -1];  % C2 (1,-sqrt(3),0)
  case 8 % #168-176, hexagonal, C6 (C6, C3h, C6h)
    % Axis convention: zC along six-fold axis, xC and yC arbitrary
    R{1} = E;
    R{2} = C2z;
    R{3} = C3zp;
    R{4} = C3zm;
    R{5} = [+1/2 -sqrt(3)/2 0; +sqrt(3)/2 +1/2 0; 0 0 +1];  % C6+ z
    R{6} = [+1/2 +sqrt(3)/2 0; -sqrt(3)/2 +1/2 0; 0 0 +1];  % C6- z
  case 9 % #177-194, hexagonal, D6 (D6, C6v, D3h, D6h)
    % Axis convention: zC along six-fold axis, xC along one of the two-fold axes
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
  case 10 % #195-206, cubic, T (T, Th)
    % Axis convention: xC, yC and zC along the three two-fold axes
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
  case 11 % #207-230, cubic, O (O, Td, Oh)
    % Axis convention: xC, yC and zC along the three four-fold axes
    R{1} = E;
    R{2} = C2z;
    R{3} = C4zp;
    R{4} = C4zm;
    R{5} = C2x;
    R{6} = C2y;
    R{7} = C2xy1;
    R{8} = C2xy2;
    R{9} = [0 0 +1; +1 0 0; 0 +1 0];  % C3+ (+1,+1,+1)
    R{10}= [0 +1 0; 0 0 +1; +1 0 0];  % C3- (+1,+1,+1)
    R{11}= [0 0 -1; -1 0 0; 0 +1 0];  % C3+ (+1,-1,-1)
    R{12}= [0 -1 0; 0 0 +1; -1 0 0];  % C3- (+1,-1,-1)
    R{13}= [0 0 +1; -1 0 0; 0 -1 0];  % C3+ (-1,+1,-1)
    R{14}= [0 -1 0; 0 0 -1; +1 0 0];  % C3- (-1,+1,-1)
    R{15}= [0 0 -1; +1 0 0; 0 -1 0];  % C3+ (-1,-1,+1)
    R{16}= [0 +1 0; 0 0 -1; -1 0 0];  % C3- (-1,-1,+1)
    R{17}= [+1 0 0; 0 0 -1; 0 +1 0];  % C4+ xC
    R{18}= [+1 0 0; 0 0 +1; 0 -1 0];  % C4- xC
    R{19}= [0 0 +1; 0 +1 0; -1 0 0];  % C4+ yC
    R{20}= [0 0 -1; 0 +1 0; +1 0 0];  % C4- yC
    R{21}= [0 0 +1; 0 -1 0; +1 0 0];  % C2 (1,0,+1)
    R{22}= [0 0 -1; 0 -1 0; -1 0 0];  % C2 (1,0,-1)
    R{23}= [-1 0 0; 0 0 +1; 0 +1 0];  % C2 (0,1,+1)
    R{24}= [-1 0 0; 0 0 -1; 0 -1 0];  % C2 (0,1,-1)
end

if transformVector
  % Apply site transformations to input vector and
  % return an array of transformed vectors
  for iR = 1:numel(R)
    vecrot(:,iR) = R{iR}*vec(:);  %#ok
  end
  out = vecrot;
else
  % Return the set of rotation matrices
  out = R;
end

end

% Determine Laue class (1-11) from space group number (1-230)
function LaueClass = groupnumber2laueclass(GroupNumber)
if     GroupNumber<=  2, LaueClass =  1; % (C1, Ci=S2), triclinic
elseif GroupNumber<= 15, LaueClass =  2; % (C2, Cs=C1h, C2h), monoclinic
elseif GroupNumber<= 74, LaueClass =  3; % (D2, C2v, D2h), orthorhombic
elseif GroupNumber<= 88, LaueClass =  4; % (C4, S4, C4h), tetragonal
elseif GroupNumber<=142, LaueClass =  5; % (D4, C4v, D2d, D4h), tetragonal
elseif GroupNumber<=148, LaueClass =  6; % (C3, C3i=S6), trigonal
elseif GroupNumber<=167, LaueClass =  7; % (D3, C3v, D3d), trigonal
elseif GroupNumber<=176, LaueClass =  8; % (C6, C3h, C6h), hexagonal
elseif GroupNumber<=194, LaueClass =  9; % (D6, C6v, D3h, D6h), hexagonal
elseif GroupNumber<=206, LaueClass = 10; % (T, Th), cubic
elseif GroupNumber<=230, LaueClass = 11; % (O, Td, Oh), cubic
end
end
