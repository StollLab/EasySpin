% symm  Determine spin Hamiltonian symmetry 
%
%   PGroup = symm(Sys)
%   [PGroup,R] = symm(Sys)
%
%   Determines the point group of the Hamiltonian
%   of a spin sytem together with its symmetry frame.
%
%   Input:
%   - Sys: Spin system specification structure
%
%   Output:
%   - PGroup: Schoenfliess point group symbol, one of
%     'Ci','C2h','D2h','C4h','D4h','S6','D3d','C6h','D6h','Th','Oh',Dinfh','O3'.
%   - R: Rotation matrix containing the axes of the
%     symmetry frame along columns.

function [PGroup,RMatrix] = symm(Sys,varargin)

if nargin==0, help(mfilename); return; end

if nargin>1
  options = varargin{end};
  DebugMode = strfind(options,'debug');
else
  DebugMode = 0;
end

[Sys,err] = validatespinsys(Sys);
error(err);
sysfields = fieldnames(Sys);


highOrderTermsPresent = ~isempty(Sys.B);

higherZeemanPresent = false;
higherzeemanFields = strncmp(sysfields,'Ham',3).';
if any(higherzeemanFields) 
  for n = find(higherzeemanFields)
    if any(Sys.(sysfields{n})(:))
        highOrderTermsPresent = true;
        higherZeemanPresent = true;
    end
  end
end

CrystalFieldPresent = false;
cf = strncmp(sysfields,'CF',2).';
if any(cf)
  for n = find(cf)
    if any(Sys.(sysfields{n})(:))
      CrystalFieldPresent = true;
      highOrderTermsPresent = true;
    end
  end
end

if DebugMode
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
  if DebugMode
    fprintf('Check whether isotropic...\n')
  end  
  if isisotropic(Sys)
    PGroup = 'O3';
    RMatrix = eye(3);
    return;
  end
end

% :TODO:
EquivalentSpins = 0;

doQMAnalysis = highOrderTermsPresent | fullTensorsGiven | EquivalentSpins;

% Geometrical analysis is flawed: doesn't work for
% CF3 radical. This has molecular symmetry C3v, should give
% spin Hamiltonian symmetry C3v x i = D3d. Geometrical
% analysis returns Ci. QM analysis gives correct D3d.

% QM analysis has problems with two axial tensors at arbitrary
% angle: returns Ci, since doesn't find common rotation axis.

if DebugMode
  if doQMAnalysis
    fprintf('Quantum mechanical symmetry analysis...\n'); 
  else
    fprintf('Geometric symmetry analysis...\n');
  end
end

if doQMAnalysis
  [PGroup, RMatrix] = symm_full(Sys,higherZeemanPresent,DebugMode);
else
  [PGroup, RMatrix] = symm_geom(Sys,DebugMode);
end

return
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Group,RMatrix] = symm_full(Sys,HigherZeemanPresent, DebugMode)

%--------------------------------------------------
% Collect all frame orientations described
% by Euler angle set from the spin system.
%--------------------------------------------------
% As defined for the spin system, AFrame and the
% like specify the Euler angles for the passive
% rotation R = erot(Sys.AFrame) which transforms
% a quantity (vector, tensor) from the molecular frame
% to the A frame. So if v is vector in
% the A frame representation, R.'*v is the same
% vector in the molecular frame representation.

% Add successively g frame, ee frame, D frame, A and Q frame
% Euler angles.
pa = [0 0 0]; % the (first) g frame itself
if isfield(Sys,'gFrame'), pa = [pa; Sys.gFrame]; end
if isfield(Sys,'eeFrame'), pa = [pa; Sys.eeFrame]; end
if isfield(Sys,'DFrame'), pa = [pa; Sys.DFrame]; end
if isfield(Sys,'AFrame')
  % make sure it works for more than 1 electron spin
  for k = 0:3:size(Sys.AFrame,2)-1
    pa = [pa; Sys.AFrame(:,k+(1:3))];
  end
end
if isfield(Sys,'QFrame'), pa = [pa; Sys.QFrame]; end

% Remove duplicates. Avoid sorting
[~,ii] = unique(pa,'rows');
pa = pa(sort(ii),:);
%--------------------------------------------------


%--------------------------------------------------
% Compute rotation matrices and z directions.
%--------------------------------------------------
% One set of Euler angles specifies exactly one
% frame, but since principal axes can be sorted
% in three different ways, each Euler angle set
% gives three potential symmetry axes and three
% symmetry frames. We generate these three by
% axes exchange (xyz) -> (yzx) -> (zxy) using
% three rotation matrices (Rz, Rx, Ry).

LockedFrame = 0;
if (~LockedFrame)
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

  % Invert frames so that all z axes point to the
  % upper hemisphere.
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
%--------------------------------------------------



%--------------------------------------------------
% Precalculation of Hamiltonian and symmetry ops.
%--------------------------------------------------
% Orientation-independent spin Hamiltonian components.
if ~HigherZeemanPresent, [F,Gx,Gy,Gz] = sham(Sys); end

% Set of field vectors for the test operations.
q = 7; % small theta aliquod
qq = 6.3456; % small phi aliquod
the = pi/180* [q,  q, q,q/2,q/2,  q,180-q,180-q,180-q,90-q,    q,  q];
phi = pi/180*([0,120,90, 90,  0,-90,    0,  -90,  180,  90,-2*qq,180]+qq);
Field = 350; % mT
FieldVecs = Field*ang2vec(phi,the); % 3x8 array
%--------------------------------------------------

if DebugMode
  fprintf('===========================================================\n');
end


%--------------------------------------------------
% Determination of point group.
%--------------------------------------------------
% At the outset we assume Ci symmetry around z axis
% of g (standard) frame. Then we try to find a higher
% symmetry in one of the frames.

iGroup = 0; % nothing
RMatrix = Rz; % eye(3)
GroupNames = {'Ci','C2h','D2h','C4h','D4h',...
  'S6','D3d','C6h','D6h','Th','Oh','Dinfh','O3'};

for iFrame = 1:nFrames % loop over all potential frames

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
  
  pg = 0;
  if HigherZeemanPresent
    eA = eig(sham(Sys,B(:,1)));
    eB = eig(sham(Sys,B(:,2)));
    eC = eig(sham(Sys,B(:,3)));
  else
    eA = eig(F + B(1,1)*Gx + B(2,1)*Gy + B(3,1)*Gz);
    eB = eig(F + B(1,2)*Gx + B(2,2)*Gy + B(3,2)*Gz);
    eC = eig(F + B(1,3)*Gx + B(2,3)*Gy + B(3,3)*Gz);

  end
  C3 = eqeig(eB,eA); % is there a C3 along z?
  C4 = eqeig(eC,eA); % is there a C4 along z?
  
  if DebugMode
    fprintf('%d C4 axes and %d C3 axes along z\n',C4,C3);
  end

  switch C4*2+C3
  case 0 % none: Ci, C2h, D2h
    if HigherZeemanPresent
      C2z = eqeig(eA,eig(sham(Sys,B(:,12))));
    else
      C2z = eqeig(eA,eig(F+B(1,12)*Gx+B(2,12)*Gy+B(3,12)*Gz));
    end
    if ~C2z
      pg = 1; % Ci
    else % D2h, C2h
      if HigherZeemanPresent, sigmaxz = eqeig(eA,eig(sham(Sys,B(:,11))));
      else sigmaxz = eqeig(eA,eig(F+B(1,11)*Gx+B(2,11)*Gy+B(3,11)*Gz));
      end
      if sigmaxz
        pg = 3; % D2h
      else
        pg = 2; % C2h
      end
    end
    
  case 1 % C3 axis: S6,D3d,Th,C6h,D6h
    if HigherZeemanPresent, sigmaxy = eqeig(eA,eig(sham(Sys,B(:,7))));
    else sigmaxy = eqeig(eA,eig(F+B(1,7)*Gx+B(2,7)*Gy+B(3,7)*Gz));
    end
    if sigmaxy % Th, C6h, D6h
      if HigherZeemanPresent,C2z = eqeig(eC,eig(sham(Sys,B(:,6))));
      else C2z = eqeig(eC,eig(F+B(1,6)*Gx+B(2,6)*Gy+B(3,6)*Gz));
      end
      if C2z % C6h, D6h
        if HigherZeemanPresent,C2x = eqeig(eC,eig(sham(Sys,B(:,8)))); 
        else C2x = eqeig(eC,eig(F+B(1,8)*Gx+B(2,8)*Gy+B(3,8)*Gz));  
        end
        if C2x, pg = 9; else pg = 8; end
      else
        pg = 10;
      end
    else % S6, D3d
      if HigherZeemanPresent,C2x = eqeig(eC,eig(sham(Sys,B(:,8))));
      else C2x = eqeig(eC,eig(F+B(1,8)*Gx+B(2,8)*Gy+B(3,8)*Gz));
      end
      if ~C2x
        if HigherZeemanPresent, C2y = eqeig(eA,eig(sham(Sys,B(:,9))));
        else C2y = eqeig(eA,eig(F+B(1,9)*Gx+B(2,9)*Gy+B(3,9)*Gz));
        end
      end
      if C2x||C2y; pg = 7; else pg = 6; end
    end
    
  case 2 % C4 axis: C4h,D4h,Oh
    
    if HigherZeemanPresent, C2x = eqeig(eC,eig(sham(Sys,B(:,8))));
    else C2x = eqeig(eC,eig(F+B(1,8)*Gx+B(2,8)*Gy+B(3,8)*Gz));
    end
    if C2x % D4h, Oh
      if HigherZeemanPresent
        Bs  = [B(2,10),B(3,10),B(1,10)];
        C3d = eqeig(eig(sham(Sys,B(:,10))),eig(sham(Sys,Bs)));
      else
        C3d = eqeig(eig(F+B(2,10)*Gx+B(3,10)*Gy+B(1,10)*Gz),...
                  eig(F+B(1,10)*Gx+B(2,10)*Gy+B(3,10)*Gz));
      end
      if C3d, pg = 11; else, pg = 5; end
    else
      pg = 4;
    end
    
  case 3 % C3 and C4 axes: Dinfh,O3
    if HigherZeemanPresent
      Cinfx = eqeig(eC,eig(sham(Sys,B(:,4))));
    else
      Cinfx = eqeig(eC,eig(F+B(1,4)*Gx+B(2,4)*Gy+B(3,4)*Gz));
    end
    if Cinfx, pg=13; else, pg=12; end
    
  end  % switch
  
  if DebugMode
    if pg>0
      fprintf('\nFrame %d: Point group %s\n',iFrame,GroupNames{pg});
    else
      fprintf('\nFrame %d: No point group identified.\n',iFrame);
    end
    fprintf('-----------------------------------------------------------\n');
  end

  % Update if symmetry is higher than current best.
  if pg>iGroup
    iGroup = pg;
    RMatrix = R;
    if DebugMode
      disp('Highest symmetry up to now!');
    end
  end
  
  if DebugMode
    disp('Columns of R: tensor axis vectors in molecular reference frame.');
    disp(R)
    if isfield(Sys,'A')
      for kk=1:size(Sys.A,1)
        if ~isfield(Sys,'AFrame')
          A = R.'*diag(Sys.A(kk,:))*R;
        else
          A = R.'*erot(Sys.AFrame(kk,:)).'*diag(Sys.A(kk,:))*erot(Sys.AFrame(kk,:))*R;
        end
        fprintf('diag(A(%d,:)) in this frame [%3.4f %3.4f %3.4f]\n',kk,A(1), ...
          A(5),A(9));
      end
    end
    pause
  end

end
%--------------------------------------------------


%--------------------------------------------------
% Output assignment.
%--------------------------------------------------
% Assertion that a point group has been found.
if iGroup==0
  error('No point group found! Save input spin system and report bug!');
end

Group = GroupNames{iGroup};
%--------------------------------------------------

if DebugMode
  fprintf('===========================================================\n');
  fprintf('Symmetry = %s\n',Group);
  %disp('A in this frame');
  %A = RMatrix.'*erot(Sys.AFrame).'*diag(Sys.A)*erot(Sys.AFrame)*RMatrix;
  %disp(A)
end

return
%-------------------------------
function really = eqeig(eA,eB)
Threshold = 1e-8;
eA = sort(eA);
eB = sort(eB);
really = norm(eA-eB) < Threshold*norm(eA);
return
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



function [SymGrp,SymFrame] = symm_geom(Sys,DebugMode)

if DebugMode
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

if DebugMode
  fprintf('  %d tensor(s)\n',numel(Sym));
end

% Remove isotropic tensors
isotropic = Sym==0;
Sym(isotropic) = [];
Ax(isotropic) = [];
Name(isotropic) = [];

if DebugMode
  fprintf('  %d anisotropic tensor(s)\n',numel(Sym));
end


% Compute total symmetry
%-----------------------------------------------------------------------
Groups = {'O3','Dinfh','D2h','C2h','Ci'};
Grp = 0;
SymFrame = eye(3);
if DebugMode, fprintf('  O3 as starting symmetry\n'); end
for iTens = 1:numel(Sym)
  [Grp,SymFrame] = combinesymms(Grp,SymFrame,Sym(iTens),Ax{iTens});
  if DebugMode
    fprintf('   + %s (%s) = %s\n',Groups{Sym(iTens)+1},Name{iTens},Groups{Grp+1});
  end
end
SymGrp = Groups{Grp+1};

return

%-------------------------------------------------------------------------------
function [SymGroup,SymFrame] = tensorsymmetry(PrincipalValues,EulerAngles)

O3 = 0; Dinfh = 1; D2h = 2;

nValueEqualities = sum(diff(sort(PrincipalValues))==0);

switch nValueEqualities
case 0
  SymGroup = D2h;
  zAxis = 3;
case 2
  SymGroup = O3;
  zAxis = 3;
case 1
  SymGroup = Dinfh;
  if (PrincipalValues(2)==PrincipalValues(3))
    zAxis = 1;
  elseif (PrincipalValues(1)==PrincipalValues(3))
    zAxis = 2;
  else 
    zAxis = 3;
  end
end

idx = [2 3 1 2 3];
SymFrame = erot(EulerAngles).';
% Columns: tensor frame axes in reference frame representation.
SymFrame = SymFrame(:,idx(zAxis+(0:2)));

return


%==========================================================================
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

function [TotalSym,TotalRot] = combinesymms(Sym1,Rot1,Sym2,Rot2)

% Rotation matrix R = erot([alpha beta gamma])
%   Cols: tensor frame axes in reference frame representation
%   Rows: reference frame axes in tensor frame representation

% Symmetry group code abbreviations
O3 = 0; Dinfh = 1; D2h = 2; C2h = 3; Ci = 4;

% Angle limits for parallel and perpendicular tests
delta = 0.2;
ParallelLimit = delta*pi/180; % maximum
PerpLimit = (90-delta)*pi/180; % minimum

% Sort Sym1 and Sym2 so that Sym1 > Sym2
% (Sym1 is of lower symmetry than Sym2)
if (Sym1<Sym2)
  [Sym2,Sym1] = deal(Sym1,Sym2);
  [Rot2,Rot1] = deal(Rot1,Rot2);
end

% any + O3, Ci + any, O3 + any
%-------------------------------------------------------------------
% - If Sym1 is Ci, the total symmetry remains Ci, since it cannot be lower.
% - If Sym2 is O3, the total symmetry is unchanged.
if (Sym1==Ci) || (Sym2==O3)
  TotalSym = Sym1;
  TotalRot = Rot1;
  return
end
% - If Sym1 is O3, the total symmetry is the new one (either Dinfh or D2h).
if (Sym1==O3)
  TotalSym = Sym2;
  TotalRot = Rot2;
  return
end

% 5 possible combinations remain:
% [1] Dinfh + Dinfh
% [2] D2h   + Dinfh
% [3] C2h   + Dinfh
% [4] D2h   + D2h
% [5] C2h   + C2h

% [1] Dinfh + Dinfh
%------------------------------------------------------------------
% a) z axes coincide -> Dinfh as before
% b) z axes perpendicular -> D2h with z perpendicular to both z1 and z2
% c) otherwise -> C2h with z perpendicular to both z1 and z2
if (Sym1==Dinfh) && (Sym2==Dinfh)
  
  z1 = Rot1(:,3);
  z2 = Rot2(:,3);
  z1z2Angle = real(acos(abs(z1'*z2)));
  if (z1z2Angle<ParallelLimit) % z axes parallel
    TotalSym = Dinfh;
    TotalRot = Rot1;
  elseif (z1z2Angle>PerpLimit) % z axes perpendicular
    TotalSym = D2h;
    newz = cross(z1,z2);
    newz = newz/norm(newz);
    TotalRot = [cross(z1,newz) z1 newz];
  else % all other angles between z1 and z2
    TotalSym = C2h;
    newz = cross(z1,z2);
    newz = newz/norm(newz);
    TotalRot = [cross(z1,newz) z1 newz];
  end

% [2] D2h + Dinfh
%------------------------------------------------------------------
% - Dinfh z parallel to any D2h axis -> D2h as before
% - Dinfh z perpendicular to one D2h axis (i.e. Dinfh z in a D2h
%   mirror plane) -> C2h with the perpendicular D2h axis as z axis
% - otherwise Ci with arbitrary orientation
elseif (Sym1==D2h) && (Sym2==Dinfh)
  z2 = Rot2(:,3); % Dinfh z axis
  z2Angles = real(acos(abs(z2'*Rot1)));
  if any(z2Angles<ParallelLimit)
    TotalSym = D2h;
    TotalRot = Rot1;
  elseif any(z2Angles>PerpLimit) % z in a D2h sigma plane
    TotalSym = C2h;
    [~,newzidx] = max(z2Angles);
    idx = [2 3 1 2 3];
    TotalRot = Rot1(:,idx(newzidx+(0:2)));
  else
    TotalSym = Ci;
    TotalRot = Rot1; % can be arbitrary
  end
  
% [3] C2h + Dinfh
%------------------------------------------------------------------
% - Dinfh z parallel to C2h z axis -> C2h as before
% - Dinfh z in the C2h mirror plane -> C2h as before
% - Dinfh z arbitrary -> Ci, arbitrary frame
elseif (Sym1==C2h) && (Sym2==Dinfh)
  
  z1 = Rot1(:,3);
  z2 = Rot2(:,3);
  z1z2Angle = real(acos(abs(z1'*z2)));
  if (z1z2Angle<ParallelLimit)
    TotalSym = C2h;
    TotalRot = Rot1;
  elseif (z1z2Angle>PerpLimit)
    TotalSym = C2h;
    TotalRot = Rot1;
  else
    TotalSym = Ci;
    TotalRot = Rot1; % can be arbitrary
  end
  
% [4] D2h + D2h
%------------------------------------------------------------------
% - all axes collinear -> D2h as before
% - only one axis collinear -> C2h with collinar axis as z
% - all axes skew -> Ci, arbitrary frame
elseif (Sym1==D2h) && (Sym2==D2h)
  
  allAngles = real(acos(abs(Rot2'*Rot1)));
  % row 1: angles of x2, row 2: angles of y2, row 3: angles of z2 
  % col 1: angles of x1, col 2: angles of y1, col 3: angles of z1
  AxisCollinear = min(allAngles)<ParallelLimit;
  if all(AxisCollinear)
    TotalSym = D2h;
    TotalRot = Rot1;
  elseif any(AxisCollinear)
    TotalSym = C2h;
    % find collinear axes
    %[minAngles,minidx2] = min(allAngles);
    [minAngles] = min(allAngles);
    [~,newzidx1] = min(minAngles);
    %newzidx2 = minidx2(newzidx1);
    % Axes newzidx2 of frame 2 and newzidx1 of frame 1 are collinear
    idx = [2 3 1 2 3];
    TotalRot = Rot1(:,idx(newzidx1+(0:2)));
  else
    TotalSym = Ci;
    TotalRot = Rot1;% can be arbitrary
  end
  
% [5] C2h + D2h
%------------------------------------------------------------------
% - any D2h axis collinar with C2h z axis -> C2h as before
% - otherwise Ci, arbitrary frame
elseif (Sym1==C2h) && (Sym2==D2h)

  z1 = Rot1(:,3);
  z1Angles = real(acos(abs(z1'*Rot2)));
  if any(z1Angles<ParallelLimit)
    TotalSym = C2h;
    TotalRot = Rot1;
  else
    TotalSym = Ci;
    TotalRot = Rot1; % can be arbitrary
  end

else

  error('Cannot compute symmetry combination: %d %d',Sym1,Sym2);

end % Combination case switchyard

return

%===============================================================================
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

return
