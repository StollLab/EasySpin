% hamsymm  Determine spin Hamiltonian symmetry 
%
%   PGroup = hamsymm(Sys)
%   [PGroup,R] = hamsymm(Sys)
%
%   Determines the point group of the Hamiltonian
%   of a spin sytem and the associated symmetry frame.
%
%   Input:
%   - Sys: Spin system structure
%
%   Output:
%   - PGroup: Schoenfliess point group symbol, one of
%     'Ci','C2h','D2h','C4h','D4h','S6','D3d','C6h','D6h','Th','Oh',Dinfh','O3'.
%   - R: Rotation matrix with the axes of the symmetry frame along columns.

function [PGroup,RMatrix] = hamsymm(Sys,varargin)

if nargin==0, help(mfilename); return; end

if nargin>1
  options = varargin{end};
  debugMode = strfind(options,'debug');
else
  debugMode = false;
end

[Sys,err] = validatespinsys(Sys);
error(err);
sysfields = fieldnames(Sys);

higherZeemanPresent = false;
higherZeemanFields = strncmp(sysfields,'Ham',3).';
if any(higherZeemanFields) 
  for n = find(higherZeemanFields)
    if any(Sys.(sysfields{n})(:))
      higherZeemanPresent = true;
    end
  end
end

crystalFieldTermsPresent = false;
cf = strncmp(sysfields,'CF',2).';
if any(cf)
  for n = find(cf)
    if any(Sys.(sysfields{n})(:))
      crystalFieldTermsPresent = true;
    end
  end
end

highOrderTermsPresent = ~isempty(Sys.B) || higherZeemanPresent || crystalFieldTermsPresent;

if debugMode
  if highOrderTermsPresent
    fprintf('High-order terms present!\n');
  else
    fprintf('No high-order terms present!\n');
  end
end

fullg = Sys.fullg;
fullA = Sys.fullA;
fullQ = Sys.fullQ;
fullD = Sys.fullD;
fullee = Sys.fullee;
fullnn = Sys.fullnn;
fullsigma = Sys.fullsigma;

fullTensorsGiven = any([fullg fullA fullD fullee fullQ fullsigma fullnn]);

if ~fullTensorsGiven && ~highOrderTermsPresent
  if debugMode
    fprintf('Check whether isotropic...\n')
  end  
  if isisotropic(Sys)
    PGroup = 'O3';
    RMatrix = eye(3);
    return;
  end
end

% :TODO:
equivalentSpins = false;

doQMAnalysis = highOrderTermsPresent || fullTensorsGiven || equivalentSpins;

% Geometrical analysis is flawed: doesn't work for CF3 radical. This has
% molecular symmetry C3v, should give spin Hamiltonian symmetry C3v x i = D3d.
% Geometrical analysis returns Ci. QM analysis gives correct D3d.

% QM analysis has problems with two axial tensors at arbitrary angle: returns
% Ci, since it is not able to identify common rotation axis.

if debugMode
  if doQMAnalysis
    fprintf('Quantum mechanical symmetry analysis...\n'); 
  else
    fprintf('Geometric symmetry analysis...\n');
  end
end

if doQMAnalysis
  [PGroup, RMatrix] = hamsymm_eigs(Sys,higherZeemanPresent,debugMode);
else
  [PGroup, RMatrix] = hamsymm_geom(Sys,debugMode);
end

% For calculations using photoselection, accurate results are in general
% only obtained if considering the full sphere (for specific combinations
% of tdm and laser polarization, lower symmetry should give an accurate
% result too, but this is safer in general).
tdmPresent = ~isempty(Sys.tdm);
if tdmPresent
  PGroup = 'C1';
end  

end
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


%-------------------------------------------------------------------------------
function [Group,RMatrix] = hamsymm_eigs(Sys,higherZeemanPresent,debugMode)

% Collect all frame orientations described by Euler angles from the spin system
%-------------------------------------------------------------------------------
% As defined for the spin system, AFrame and the like specify the Euler angles
% for the passive rotation R_M2A = erot(Sys.AFrame) which transforms a quantity
% (vector, tensor) from the molecular frame to the A frame. To go from the A
% frame to the molecular frame, use R_A2M = R_M2A.'. If v_A is vector in
% the A frame representation, v_M = R_M2A.'*v_A is the same vector in the
% molecular frame representation. For tensors: T_M = R_M2A.'*T_A*R_M2A.

% Add successively g frame, ee frame, D frame, A and Q frame Euler angles.
pa = [0 0 0];  % the molecular frame
if isfield(Sys,'gFrame'), pa = [pa; Sys.gFrame]; end
if isfield(Sys,'eeFrame'), pa = [pa; Sys.eeFrame]; end
if isfield(Sys,'DFrame'), pa = [pa; Sys.DFrame]; end
if isfield(Sys,'AFrame')
  % make sure it works for more than 1 electron spin
  pa = [pa; reshape(Sys.AFrame.',3,[]).'];
end
if isfield(Sys,'QFrame'), pa = [pa; Sys.QFrame]; end
if isfield(Sys,'nnFrame'), pa = [pa; Sys.nnFrame]; end
if isfield(Sys,'sigmaFrame'), pa = [pa; Sys.sigmaFrame]; end

% Remove duplicates and sort
[~,ii] = unique(pa,'rows');
pa = pa(sort(ii),:);


% Compute rotation matrices and z directions
%-------------------------------------------------------------------------------
% One set of Euler angles specifies exactly one frame, but since principal axes
% can be sorted in three different ways, each Euler angle set gives three
% potential symmetry axes and three symmetry frames. We generate these three by
% axes exchange (xyz) -> (yzx) -> (zxy) using three rotation matrices
% (Rz, Rx, Ry).

lockedFrame = 0;
if ~lockedFrame
  % Initialize.
  Rz = eye(3);
  Rx = Rz([2 3 1],:);
  Ry = Rz([3 1 2],:);

  nFrames = 3*size(pa,1);

  % Rotation matrices.
  Rots = zeros(3,9*size(pa,1));
  for k = 1:size(pa,1)
    R = erot(pa(k,:)).';
    Rots(:,k*9+(-8:0)) = [R*Rz, R*Rx, R*Ry];
  end
  Rots = reshape(Rots,3,3,nFrames);

  % Invert frames so that all z axes point to the upper hemisphere.
  %Invert = Rots(3,3,:)<0;
  %Rots(:,:,Invert) = -Rots(:,:,Invert);

  % Remove duplicate frames. Avoid sorting.
  [~,idx] = unique(reshape(Rots,9,nFrames).','rows');
  Rots = Rots(:,:,sort(idx));
  nFrames = size(Rots,3);
else
  Rots = eye(3);
  nFrames = 1;
end



% Precalculation of Hamiltonian and symmetry ops
%-------------------------------------------------------------------------------
% Orientation-independent spin Hamiltonian components.
if ~higherZeemanPresent
  [H0,mux,muy,muz] = ham(Sys);
end

% Set of field vectors for the test operations.
q = 7.3234;  % small theta increment
qq = 6.3456;  % small phi increment
the = pi/180* [q,  q, q,q/2,q/2,  q,180-q,180-q,180-q,90-q,    q,  q];
phi = pi/180*([0,120,90, 90,  0,-90,    0,  -90,  180,  90,-2*qq,180]+qq);
Field = 350; % mT
FieldVecs = Field*ang2vec(phi,the);  % 3x12 array

if debugMode
  fprintf('===========================================================\n');
end


% Determination of point group
%-------------------------------------------------------------------------------
% At the outset we assume Ci symmetry around z axis
% of g (standard) frame. Then we try to find a higher
% symmetry in one of the frames.

highestSymmetry = 0;  % no symmetry (i.e. C1)
RMatrix = Rz;  % eye(3)
GroupNames = {'Ci','C2h','D2h','C4h','D4h',...
  'S6','D3d','C6h','D6h','Th','Oh','Dinfh','O3'};

if debugMode
  fprintf('initial symmetry: C1\n');
end

for iFrame = 1:nFrames % loop over all potential frames
  if debugMode
    fprintf('Frame %d ------------------------\n',iFrame)
  end

  % Rotate all orientations.
  % This is a passive rotation of the FieldVecs from
  % the current frame to the standard frame, or, what
  % is the point here, the ACTIVE rotation of FieldVecs:
  % A vector v = [0;0;1] (understood in the standard frame)
  % will point along the z axis of the current frame after
  % the transformation v1 = R*v. The vector's orientation
  % has changed, but not its representation.
  R = Rots(:,:,iFrame);
  B = R*FieldVecs;

  if higherZeemanPresent
    eA = eig(ham(Sys,B(:,1)));
    eB = eig(ham(Sys,B(:,2)));
    eC = eig(ham(Sys,B(:,3)));
  else
    eA = eig(H0 - B(1,1)*mux - B(2,1)*muy - B(3,1)*muz);
    eB = eig(H0 - B(1,2)*mux - B(2,2)*muy - B(3,2)*muz);
    eC = eig(H0 - B(1,3)*mux - B(2,3)*muy - B(3,3)*muz);
  end

  C3 = eqeig(eB,eA);  % is there a C3 along z?
  C4 = eqeig(eC,eA);  % is there a C4 along z?
  
  if debugMode
    if C3
      fprintf('C3z axis found\n');
    else
      fprintf('No C3z axis found\n');
    end
    if C4
      fprintf('C4z axis found\n');
    else
      fprintf('No C4z axis found\n');
    end
  end

  if ~C3 && ~C4  % none: Ci, C2h, D2h

    if higherZeemanPresent
      C2z = eqeig(eA,eig(ham(Sys,B(:,12))));
    else
      C2z = eqeig(eA,eig(H0-B(1,12)*mux-B(2,12)*muy-B(3,12)*muz));
    end
    if ~C2z
      pg = 1;  % Ci
    else % D2h, C2h
      if higherZeemanPresent
        sigmaxz = eqeig(eA,eig(ham(Sys,B(:,11))));
      else
        sigmaxz = eqeig(eA,eig(H0-B(1,11)*mux-B(2,11)*muy-B(3,11)*muz));
      end
      if sigmaxz
        pg = 3;  % D2h
      else
        pg = 2;  % C2h
      end
    end
    
  elseif C3 && ~C4  % C3 axis: S6, D3d, Th, C6h, D6h

    if higherZeemanPresent
      sigmaxy = eqeig(eA,eig(ham(Sys,B(:,7))));
    else
      sigmaxy = eqeig(eA,eig(H0-B(1,7)*mux-B(2,7)*muy-B(3,7)*muz));
    end
    if sigmaxy  % Th, C6h, D6h
      if higherZeemanPresent
        C2z = eqeig(eC,eig(ham(Sys,B(:,6))));
      else
        C2z = eqeig(eC,eig(H0-B(1,6)*mux-B(2,6)*muy-B(3,6)*muz));
      end
      if C2z  % C6h, D6h
        if higherZeemanPresent
          C2x = eqeig(eC,eig(ham(Sys,B(:,8)))); 
        else
          C2x = eqeig(eC,eig(H0-B(1,8)*mux-B(2,8)*muy-B(3,8)*muz));  
        end
        if C2x
          pg = 9;  % D6h
        else
          pg = 8;  % C6h
        end
      else
        pg = 10;  % Th
      end
    else  % S6, D3d
      if higherZeemanPresent
        C2x = eqeig(eC,eig(ham(Sys,B(:,8))));
      else
        C2x = eqeig(eC,eig(H0-B(1,8)*mux-B(2,8)*muy-B(3,8)*muz));
      end
      if ~C2x
        if higherZeemanPresent
          C2y = eqeig(eA,eig(ham(Sys,B(:,9))));
        else
          C2y = eqeig(eA,eig(H0-B(1,9)*mux-B(2,9)*muy-B(3,9)*muz));
        end
      end
      if C2x || C2y
        pg = 7;  % D3d
      else
        pg = 6;  % S6
      end
    end
    
  elseif ~C3 && C4  % C4 axis: C4h, D4h, Oh
    
    if higherZeemanPresent
      C2x = eqeig(eC,eig(ham(Sys,B(:,8))));
    else
      C2x = eqeig(eC,eig(H0-B(1,8)*mux-B(2,8)*muy-B(3,8)*muz));
    end
    if C2x % D4h, Oh
      if higherZeemanPresent
        Bs  = [B(2,10),B(3,10),B(1,10)];
        C3d = eqeig(eig(ham(Sys,B(:,10))),eig(ham(Sys,Bs)));
      else
        C3d = eqeig(eig(H0-B(2,10)*mux-B(3,10)*muy-B(1,10)*muz),...
                    eig(H0-B(1,10)*mux-B(2,10)*muy-B(3,10)*muz));
      end
      if C3d
        pg = 11;  % Oh
      else
        pg = 5;  % D4h
      end
    else
      pg = 4;  % C4h
    end
    
  else  % C3 and C4 axes: Dinfh, O3

    if higherZeemanPresent
      Cinfx = eqeig(eC,eig(ham(Sys,B(:,4))));
    else
      Cinfx = eqeig(eC,eig(H0-B(1,4)*mux-B(2,4)*muy-B(3,4)*muz));
    end
    if Cinfx
      pg = 13;  % O3
    else
      pg = 12;  % Dinfh
    end
    
  end
  
  if debugMode
    if pg>0
      fprintf('Point group %s\n',GroupNames{pg});
    else
      fprintf('No point group identified.\n');
    end
  end

  % Update if symmetry is higher than current best.
  if pg>highestSymmetry
    highestSymmetry = pg;
    RMatrix = R;
    if debugMode
      disp('Highest symmetry up to now.');
    end
  end
  
end


% Output assignment
%-------------------------------------------------------------------------------
% Assertion that a point group has been found.
if highestSymmetry==0
  error('No point group found! Save input spin system and report bug!');
end

Group = GroupNames{highestSymmetry};

if debugMode
  fprintf('===========================================================\n');
  fprintf('Symmetry = %s\n',Group);
end

end
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


%===============================================================================
% Determine if two sets of eigenvalues are (approximately) equal
function tf = eqeig(eA,eB)
threshold = 1e-12;
eA = sort(eA);
eB = sort(eB);
tf = norm(eA-eB) < threshold*norm(eA);
end


%===============================================================================
% Determine point group symmetry using geometric reasoning
function [SymGrp,symFrame] = hamsymm_geom(Sys,debugMode)

if debugMode
  fprintf('===========================================================\n');
end

Sys = validatespinsys(Sys);
nElectrons = Sys.nElectrons;
nNuclei = Sys.nNuclei;
nElCouplings = nElectrons*(nElectrons-1)/2;
nNucCouplings = nNuclei*(nNuclei-1)/2;

% Determine symmetries of all tensors in spin system
%-------------------------------------------------------------------------------
Sym = [];
Ax = {};
Name = {};

% Electron Zeeman interaction
for iE = 1:nElectrons
  [Sym(end+1),Ax{end+1}] = tensorsymmetry(Sys.g(iE,:),Sys.gFrame(iE,:));
  Name{end+1} = sprintf('g%d',iE);
end

% Zero-field interaction
for iE = 1:nElectrons
  if Sys.S(iE)>1/2
    [Sym(end+1),Ax{end+1}] = tensorsymmetry(Sys.D(iE,:),Sys.DFrame(iE,:));
    Name{end+1} = sprintf('D%d',iE);
  end
end

% Electron-electron interaction
for iC = 1:nElCouplings
  [Sym(end+1),Ax{end+1}] = tensorsymmetry(Sys.ee(iC,:),Sys.eeFrame(iC,:));
  Name{end+1} = sprintf('ee%d',iC);
end

% Nuclear Zeeman interaction
for iN = 1:nNuclei
  [Sym(end+1),Ax{end+1}] = tensorsymmetry(Sys.sigma(iN,:),Sys.sigmaFrame(iN,:));
  Name{end+1} = sprintf('sigma%d',iE);
end

% Hyperfine interaction
for iN = 1:nNuclei
  eidx = 1:3;
  for iE = 1:nElectrons
    [Sym(end+1),Ax{end+1}] = tensorsymmetry(Sys.A(iN,eidx),Sys.AFrame(iN,eidx));
    Name{end+1} = sprintf('A%d%d',iE,iN);
    eidx = eidx + 3;
  end
end

% Nuclear quadrupole interaction
for iN = 1:nNuclei
  if Sys.I(iN)>1/2
    [Sym(end+1),Ax{end+1}] = tensorsymmetry(Sys.Q(iN,:),Sys.QFrame(iN,:));
    Name{end+1} = sprintf('Q%d',iN);
  end
end

% Nucleus-nucleus interaction
for iC = 1:nNucCouplings
  [Sym(end+1),Ax{end+1}] = tensorsymmetry(Sys.nn(iC,:),Sys.nnFrame(iC,:));
  Name{end+1} = sprintf('nn%d',iC);
end

if debugMode
  fprintf('  %d tensor(s)\n',numel(Sym));
end

% Remove isotropic tensors
isotropic = Sym==0;
Sym(isotropic) = [];
Ax(isotropic) = [];
Name(isotropic) = [];

if debugMode
  fprintf('  %d anisotropic tensor(s)\n',numel(Sym));
end


% Compute total symmetry
%-------------------------------------------------------------------------------
pointGroups = {'O3','Dinfh','D2h','C2h','Ci'};
Grp = 0;
symFrame = eye(3);
if debugMode, fprintf('  O3 as starting symmetry\n'); end
for iTens = 1:numel(Sym)
  [Grp,symFrame] = combinesymms(Grp,symFrame,Sym(iTens),Ax{iTens});
  if debugMode
    fprintf('   + %s (%s) = %s\n',pointGroups{Sym(iTens)+1},Name{iTens},pointGroups{Grp+1});
  end
end
SymGrp = pointGroups{Grp+1};

end


%===============================================================================
function [symGroup,symFrame] = tensorsymmetry(principalValues,EulerAngles)

O3 = 0;
Dinfh = 1;
D2h = 2;

nValueEqualities = sum(diff(sort(principalValues))==0);

switch nValueEqualities
case 0
  symGroup = D2h;
  zAxis = 3;
case 2
  symGroup = O3;
  zAxis = 3;
case 1
  symGroup = Dinfh;
  if principalValues(2)==principalValues(3)
    zAxis = 1;
  elseif principalValues(1)==principalValues(3)
    zAxis = 2;
  else 
    zAxis = 3;
  end
end

idx = [2 3 1 2 3];
symFrame = erot(EulerAngles).';
% Columns: tensor frame axes in reference frame representation.
symFrame = symFrame(:,idx(zAxis+(0:2)));

end


%===============================================================================
% combinesymms  Symmetry group combination
%
%      [Sym,Rot] = combinesymms(Sym1,Rot1,Sym2,Rot2)
%
%      Computes the combined symmetry group of Sym1 and Sym2.
%      Sym1 and Sym2 are symmetry group IDs: 0 O3, 1 Dinfh, 2 D2h,
%      3 C2h, 4 Ci. Rot1 and Rot2 are the rotation matrices containing
%      the orientations of the principal symmetry axes in a reference frame
%      Sym is the combined symmetry group, and Rot contains the combined axes
%      orientation. Sym2 must be either O3, Dinf or D2h.

function [totalSym,totalRot] = combinesymms(Sym1,Rot1,Sym2,Rot2)

% Rotation matrix R = erot([alpha beta gamma])
%   Cols: tensor frame axes in reference frame representation
%   Rows: reference frame axes in tensor frame representation

% Symmetry group code abbreviations
O3 = 0; Dinfh = 1; D2h = 2; C2h = 3; Ci = 4;

% Angle limits for parallel and perpendicular tests
delta = 0.2;  % deg
parallelLimit = delta*pi/180;  % maximum
perpLimit = (90-delta)*pi/180;  % minimum

% Sort Sym1 and Sym2 so that Sym1 > Sym2
% (Sym1 is of lower symmetry than Sym2)
if Sym1<Sym2
  [Sym2,Sym1] = deal(Sym1,Sym2);
  [Rot2,Rot1] = deal(Rot1,Rot2);
end

% any + O3, Ci + any, O3 + any
%-------------------------------------------------------------------------------
% - If Sym1 is Ci, the total symmetry remains Ci, since it cannot be lower.
% - If Sym2 is O3, the total symmetry is unchanged.
if Sym1==Ci || Sym2==O3
  totalSym = Sym1;
  totalRot = Rot1;
  return
end
% - If Sym1 is O3, the total symmetry is the new one (either Dinfh or D2h).
if Sym1==O3
  totalSym = Sym2;
  totalRot = Rot2;
  return
end

% 5 possible combinations remain:
% [1] Dinfh + Dinfh
% [2] D2h   + Dinfh
% [3] C2h   + Dinfh
% [4] D2h   + D2h
% [5] C2h   + C2h

% [1] Dinfh + Dinfh
%-------------------------------------------------------------------------------
% a) z axes coincide -> Dinfh as before
% b) z axes perpendicular -> D2h with z perpendicular to both z1 and z2
% c) otherwise -> C2h with z perpendicular to both z1 and z2
if Sym1==Dinfh && Sym2==Dinfh
  
  z1 = Rot1(:,3);
  z2 = Rot2(:,3);
  z1z2Angle = real(acos(abs(z1'*z2)));
  if z1z2Angle<parallelLimit  % z axes parallel
    totalSym = Dinfh;
    totalRot = Rot1;
  elseif z1z2Angle>perpLimit  % z axes perpendicular
    totalSym = D2h;
    newz = cross(z1,z2);
    newz = newz/norm(newz);
    totalRot = [cross(z1,newz) z1 newz];
  else % all other angles between z1 and z2
    totalSym = C2h;
    newz = cross(z1,z2);
    newz = newz/norm(newz);
    totalRot = [cross(z1,newz) z1 newz];
  end

% [2] D2h + Dinfh
%-------------------------------------------------------------------------------
% - Dinfh z parallel to any D2h axis -> D2h as before
% - Dinfh z perpendicular to one D2h axis (i.e. Dinfh z in a D2h
%   mirror plane) -> C2h with the perpendicular D2h axis as z axis
% - otherwise Ci with arbitrary orientation
elseif Sym1==D2h && Sym2==Dinfh
  z2 = Rot2(:,3); % Dinfh z axis
  z2Angles = real(acos(abs(z2'*Rot1)));
  if any(z2Angles<parallelLimit)
    totalSym = D2h;
    totalRot = Rot1;
  elseif any(z2Angles>perpLimit) % z in a D2h sigma plane
    totalSym = C2h;
    [~,newzidx] = max(z2Angles);
    idx = [2 3 1 2 3];
    totalRot = Rot1(:,idx(newzidx+(0:2)));
  else
    totalSym = Ci;
    totalRot = Rot1; % can be arbitrary
  end
  
% [3] C2h + Dinfh
%-------------------------------------------------------------------------------
% - Dinfh z parallel to C2h z axis -> C2h as before
% - Dinfh z in the C2h mirror plane -> C2h as before
% - Dinfh z arbitrary -> Ci, arbitrary frame
elseif Sym1==C2h && Sym2==Dinfh
  
  z1 = Rot1(:,3);
  z2 = Rot2(:,3);
  z1z2Angle = real(acos(abs(z1'*z2)));
  if z1z2Angle<parallelLimit
    totalSym = C2h;
    totalRot = Rot1;
  elseif z1z2Angle>perpLimit
    totalSym = C2h;
    totalRot = Rot1;
  else
    totalSym = Ci;
    totalRot = Rot1; % can be arbitrary
  end
  
% [4] D2h + D2h
%-------------------------------------------------------------------------------
% - all axes collinear -> D2h as before
% - only one axis collinear -> C2h with collinar axis as z
% - all axes skew -> Ci, arbitrary frame
elseif Sym1==D2h && Sym2==D2h
  
  allAngles = real(acos(abs(Rot2'*Rot1)));
  % row 1: angles of x2, row 2: angles of y2, row 3: angles of z2 
  % col 1: angles of x1, col 2: angles of y1, col 3: angles of z1
  AxisCollinear = min(allAngles)<parallelLimit;
  if all(AxisCollinear)
    totalSym = D2h;
    totalRot = Rot1;
  elseif any(AxisCollinear)
    totalSym = C2h;
    % find collinear axes
    %[minAngles,minidx2] = min(allAngles);
    [minAngles] = min(allAngles);
    [~,newzidx1] = min(minAngles);
    %newzidx2 = minidx2(newzidx1);
    % Axes newzidx2 of frame 2 and newzidx1 of frame 1 are collinear
    idx = [2 3 1 2 3];
    totalRot = Rot1(:,idx(newzidx1+(0:2)));
  else
    totalSym = Ci;
    totalRot = Rot1;% can be arbitrary
  end
  
% [5] C2h + D2h
%-------------------------------------------------------------------------------
% - any D2h axis collinar with C2h z axis -> C2h as before
% - otherwise Ci, arbitrary frame
elseif Sym1==C2h && Sym2==D2h

  z1 = Rot1(:,3);
  z1Angles = real(acos(abs(z1'*Rot2)));
  if any(z1Angles<parallelLimit)
    totalSym = C2h;
    totalRot = Rot1;
  else
    totalSym = Ci;
    totalRot = Rot1; % can be arbitrary
  end

else

  error('Cannot compute symmetry combination: %d %d',Sym1,Sym2);

end % Combination case switchyard

end

%===============================================================================
% Determine if the spin system contains isotropic terms only
function iso = isisotropic(Sys)

% fn = fieldnames(Sys);
% HighOrderTerm = strncmp(fn,'B',1);
% for k = 1:numel(fn)
%   HighOrderTerm(k) = HighOrderTerm(k) & (3==numel(fn(k)));
% end
% if any(HighOrderTerm), return; end

isiso = @(T) isfield(Sys,T) && any(any(diff(Sys.(T),[],2)));

iso = isiso('g') && isiso('D') && isiso('ee') && ...
      isiso('Q') && ...
      isiso('sigma') && isiso('nn');

if iso && isfield(Sys,'A')
  eidx = 1:3;
  for e = 1:Sys.nElectrons
    iso = iso && any(diff(Sys.A(:,eidx),[],2));
    eidx = eidx + 3;
  end
end

end
