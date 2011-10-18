function out = sitetransforms(ID,x)

% ID: space group symbol or number, or point group symbol
% x:  vector(s) to transform (optional)

persistent SpaceGroups
if isempty(SpaceGroups)
  esPath = fileparts(which(mfilename));
  DataFile = [esPath filesep 'spacegroups.txt'];
  SpaceGroups = textread(DataFile,'%s');
end
PointGroupsS = {'C1','Ci','C2','Cs','C2h','D2','C2v','D2h','C4','S4','C4h',...
  'D4','C4v','D2d','D4h','C3','C3i','D3','C3v','D3d','C6','C3h','C6h',...
  'D6','C6v','D3h','D6h','T','Th','O','Td','Oh'};
PointGroupsHM = {'1','-1','2','m','2/m','222','mm2','mmm','4','-4','4/m',...
  '422','4mm','-42m','4/mmm','3','-3','32','3m','-3m','6','-6','6/m',...
  '622','6mm','-6m2','6/mmm','23','m-3','432','-43m','m-3m'};
LaueString = {'L1','L2','L3','L4','L5','L6','L7','L8','L9','L10','L11'};
Laue = [1 1 2 2 2 3 3 3 4 4 4 5 5 5 5 6 6 7 7 7 8 8 8 9 9 9 9 10 10 11 11 11];

% Determine Laue class from space group number, space group
% symbol, or point group symbol
%--------------------------------------------------------
LaueClass = [];
if ischar(ID)
  idx = strmatch(ID,PointGroupsS,'exact');
  if isempty(idx)
    idx = strmatch(ID,PointGroupsHM,'exact');
  end
  if ~isempty(idx)
    LaueClass = Laue(idx);
  else
    idx = strmatch(ID,SpaceGroups,'exact');
    if ~isempty(idx)
      ID = idx;
    else
      idx = strmatch(ID,LaueString,'exact');
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
    error('Space group number be between 1 and 230.');
  elseif mod(ID,1)~=0 || ~isreal(ID) || (ID<1) || (ID>230)
    error('Wrong space group number %g. Must be between 1 and 230.',ID);
  end
  if     ID<=  2, LaueClass =  1;
  elseif ID<= 15, LaueClass =  2;
  elseif ID<= 74, LaueClass =  3;
  elseif ID<= 88, LaueClass =  4;
  elseif ID<=142, LaueClass =  5;
  elseif ID<=148, LaueClass =  6;
  elseif ID<=167, LaueClass =  7;
  elseif ID<=176, LaueClass =  8;
  elseif ID<=194, LaueClass =  9;
  elseif ID<=206, LaueClass = 10;
  elseif ID<=230, LaueClass = 11;
  end
end

% Compile list of rotational symmetry transformation
% for the given Laue class
%-----------------------------------------------------------
E = [1 0 0, 0 1 0, 0 0 1];
C2x = [1 0 0, 0 -1 0, 0 0 -1];
C2y = [-1 0 0, 0 1 0, 0 0 -1];
C2z = [-1 0 0, 0 -1 0, 0 0 1];
C4zp = [0 -1 0, +1 0 0, 0 0 1];
C4zm = [0 +1 0, -1 0 0, 0 0 1];
C2d1 = [0 +1 0, +1 0 0, 0 0 -1]; % around (x,y,0)
C2d2 = [0 -1 0, -1 0 0, 0 0 -1]; % around (x,-y,0)
C3zp = [-1/2 -sqrt(3)/2 0, +sqrt(3)/2 -1/2 0, 0 0 1];
C3zm = [-1/2 +sqrt(3)/2 0, -sqrt(3)/2 -1/2 0, 0 0 1];

switch LaueClass
  case 1, % Ci (C1 etc)
    R = E;
  case 2, % C2h (C2 etc)
    R = [E; C2z];
  case 3,% D2h (D2 etc)
    R = [E; C2x; C2y; C2z];
  case 4, % C4h (C4 etc)
    R = [E; C2z; C4zp; C4zm];
  case 5, % D4h (D4 etc)
    R = [E; C2x; C2y; C2z; C4zp; C4zm; C2d1; C2d2];
  case 6, % C3i (C3 etc)
    R = [E; C3zp; C3zm];
  case 7, % D3d (D3 etc)
    R = [E; C3zp; C3zm; C2x;
          -1/2 +sqrt(3)/2 0, +sqrt(3)/2 +1/2 0, 0 0 -1;  % C2(x,+sqrt(3)y)
          -1/2 -sqrt(3)/2 0, -sqrt(3)/2 +1/2 0, 0 0 -1]; % C2(x,-sqrt(3)y)
  case 8, % C6h (C6 etc)
    R = [E; C3zp; C3zm; C2z;
         +1/2 -sqrt(3)/2 0, +sqrt(3)/2 +1/2 0, 0 0 1; % C6(z)+
         +1/2 +sqrt(3)/2 0, -sqrt(3)/2 +1/2 0, 0 0 1]; % C6(z)- 
  case 9, % D6h (D6 etc)
    R = [E; C3zp; C3zm; C2z; C2x; C2y;
         +1/2 -sqrt(3)/2 0, +sqrt(3)/2 +1/2 0, 0 0 1; % C6(z)+
         +1/2 +sqrt(3)/2 0, -sqrt(3)/2 +1/2 0, 0 0 1; % C6(z)- 
          -1/2 +sqrt(3)/2 0, +sqrt(3)/2 +1/2 0, 0 0 -1; % C2(x,+sqrt(3)y)
          -1/2 -sqrt(3)/2 0, -sqrt(3)/2 +1/2 0, 0 0 -1; % C2(x,-sqrt(3)y)
          +1/2 +sqrt(3)/2 0, +sqrt(3)/2 -1/2 0, 0 0 -1; % C2(sqrt(3)x,+y)
          +1/2 -sqrt(3)/2 0, -sqrt(3)/2 -1/2 0, 0 0 -1]; % C2(sqrt(3)x,-y)
  case 10, % Th (T etc)
    R = [E; C2x; C2y; C2z;
          0 0 +1, +1 0 0, 0 1 0; % C3(+++)+
          0 0 -1, -1 0 0, 0 1 0; % C3(+--)+
          0 0 +1, -1 0 0, 0 -1 0; % C3(-+-)+
          0 0 -1, +1 0 0, 0 -1 0; % C3(--+)+
          0 +1 0, 0 0 1, +1 0 0; % C3(+++)-
          0 -1 0, 0 0 1, -1 0 0; % C3(+--)-
          0 -1 0, 0 0 -1, +1 0 0; % C3(-+-)-
          0 +1 0, 0 0 -1, -1 0 0]; % C3(--+)-
  case 11, % Oh (O etc)
    R = [E; C2x; C2y; C2z; C4zp; C4zm; C2d1; C2d2;
          0 0 +1, +1 0 0, 0 1 0; % C3(+++)+
          0 0 -1, -1 0 0, 0 1 0; % C3(+--)+
          0 0 +1, -1 0 0, 0 -1 0; % C3(-+-)+
          0 0 -1, +1 0 0, 0 -1 0; % C3(--+)+
          0 +1 0, 0 0 1, +1 0 0; % C3(+++)-
          0 -1 0, 0 0 1, -1 0 0; % C3(+--)-
          0 -1 0, 0 0 -1, +1 0 0; % C3(-+-)-
          0 +1 0, 0 0 -1, -1 0 0; % C3(--+)-
          1 0 0, 0 0 -1, 0 +1 0; % C4x+
          1 0 0, 0 0 +1, 0 -1 0; % C4x-
          0 0 +1, 0 1 0, -1 0 0; % C4y+
          0 0 -1, 0 1 0, +1 0 0; % C4y-
          0 0 +1, 0 -1 0, +1 0 0; % C2 (101)
          0 0 -1, 0 -1 0, -1 0 0; % C2 (10-1)
          -1 0 0, 0 0 +1, 0 +1 0; % C2 (011)
          -1 0 0, 0 0 -1, 0 -1 0]; % C2 (01-1)
end

nR = size(R,1);

for k=1:nR
  R_{k} = reshape(R(k,:),3,3).';
end
R = R_;

if (nargin>1)
  for k=1:nR
    v(:,k) = R{k}*x(:);
  end
  out = v;
else
  out = R;
end
