% hamsymm  Determine spin Hamiltonian symmetry 
%
%   PGroup = hamsymm(Sys)
%   [PGroup,R_S2M] = hamsymm(Sys)
%
%   Determines the point group of the Hamiltonian of a spin sytem and the
%   associated symmetry frame.
%
%   Input:
%     Sys     Spin system structure
%
%   Output:
%     PGroup  Schoenfliess point group symbol, one of
%              'Ci','C2h','D2h','C4h','D4h','S6','D3d',
%              'C6h','D6h','Th','Oh','Dinfh','O3'.
%     R_S2M   Transformation matrix from the symmetry frame to the molecular frame.

% This function implements two algorithms, one based on simple geometric
% comparison and combination of second-order tensor symmetries, and one
% more general (but more expensive) one based on comparing spin Hamiltonian
% eigenvalues. Neither is perfect.
%
% The geometric analysis disregards tensor principal values and it is therefore
% not able to identify higher symmetries when multiple tensors with identical
% principal values but different orientations are present. For example, it
% doesn't work for the CF3 radical, with three identical A tensors rotated
% around the zM axis. This has molecular symmetry C3v and spin Hamiltonian
% symmetry C3v x i = D3d. Geometric analysis returns Ci.
%
%   A = [80 87 262]*2.8;
%   Sys.A = [A; A; A];
%   beta = 17.8;  % deg
%   Sys.AFrame = deg2rad([0 beta 0; -120 beta 0; 120 beta 0]);
%
% The eigenvalue analysis fails in any situation where the symmetry frame
% is different from any tensor frame, e.g. two identical axial hf tensors
% with AFrame = deg2rad([20 30 17; 180+20 30 17]).
%
% Thankfully, "failure" results in a lower symmetry than the actual one
% being identified, so the simulated spectrum is still correct.
%
% A possible route for improvement is to generate more candidate symmetry
% frames when multiple equivalent tensors are present. Also, a more general
% geometric analysis could follow the same logic as the eigenvalue
% approach, by examining the tensors as tri-axial ellispoids.

function [PGroup,RMatrix] = hamsymm(Sys)

if nargin==0, help(mfilename); return; end

[Sys,err] = validatespinsys(Sys);
error(err);

% For calculations using photoselection, accurate results are in general only
% obtained if considering the full sphere (for specific combinations of tdm and
% laser polarization, lower symmetry should give an accurate result too, but
% this is safer in general).
tdmPresent = ~isempty(Sys.tdm);
if tdmPresent
  PGroup = 'C1';
  RMatrix = eye(3);
  return
end

% Check for complicating circumstances (anything beyond standard tensors)
higherOrderZeeman = checkFields(Sys,'Ham');
crystalFieldTerms = checkFields(Sys,'CF');
StevensTerms = ~isempty(Sys.B);
fullTensors = any([Sys.fullg Sys.fullA Sys.fullD Sys.fullee Sys.fullQ  Sys.fullsigma Sys.fullnn]);
onlySimpleTensors = ~StevensTerms && ~higherOrderZeeman &&  ~crystalFieldTerms && ~fullTensors;

% Determine symmetry
if onlySimpleTensors
  if isisotropic(Sys)
    PGroup = 'O3';
    RMatrix = eye(3);
  else
    [PGroup, RMatrix] = hamsymm_geom(Sys);
  end
else
  [PGroup, RMatrix] = hamsymm_eigs(Sys,higherOrderZeeman);
end
RMatrix = RMatrix.';  % conver to molecular frame -> symmetry frame

end
%===============================================================================


%-------------------------------------------------------------------------------
function [Group,RMatrix] = hamsymm_eigs(Sys,higherZeemanPresent)

% Collect all tensor frame orientations from the spin system
%-------------------------------------------------------------------------------
% As defined for the spin system, AFrame and the like specify the Euler angles
% for the passive rotation R_M2A = erot(Sys.AFrame) which transforms a quantity
% (vector, tensor) from the molecular frame to the A frame. To go from the A
% frame to the molecular frame, use R_A2M = R_M2A.'. If v_A is vector in
% the A frame representation, v_M = R_A2M*v_A is the same vector in the
% molecular frame representation. For tensors: T_M = R_A2M*T_A*R_A2M.'.

% Add frames from g, ee, D, A, Q, nn and sigma (which all have 3 columns)
frameAngles = [0 0 0];  % the molecular frame
fnames = ["gFrame","eeFrame","DFrame","QFrame","nnFrame","sigmaFrame"];
for fn = fnames
  if isfield(Sys,fn)
    frameAngles = [frameAngles; Sys.(fn)];
  end
end
% Separate treatment for Sys.A, since it might have more than 3 columns
if isfield(Sys,'AFrame')
  frameAngles = [frameAngles; reshape(Sys.AFrame.',3,[]).'];
end

% Remove duplicates and sort
[~,idx] = unique(frameAngles,'rows');
frameAngles = frameAngles(sort(idx),:);


% Calculate rotation matrices for all symmetry frames
%-------------------------------------------------------------------------------
% One set of Euler angles specifies exactly one frame, but since principal axes
% can be sorted in three different ways, each Euler angle set gives three
% potential symmetry axes and three symmetry frames. We generate these three by
% axes exchange (xyz) -> (yzx) -> (zxy) using three rotation matrices
% (Rz, Rx, Ry).

% Rotation matrices for all symmetry frames (including axis permutations)
nFrames = size(frameAngles,1);
allFrames = zeros(3,3,3*nFrames);
I = eye(3);  % permutation 123 -> 123
R231 = I([2 3 1],:);  % permutation 123 -> 231
R312 = I([3 1 2],:);  % permutation 123 -> 312
for f = 1:nFrames
  R_M2S = erot(frameAngles(f,:));
  R_S2M = R_M2S.';
  allFrames(:,:,(f-1)*3+(1:3)) = reshape(R_S2M*[I, R231, R312],3,3,3);
end

% Remove duplicate frames
[~,idx] = unique(reshape(allFrames,9,[]).','rows');
allFrames = allFrames(:,:,sort(idx));


% Precalculation of Hamiltonian
%-------------------------------------------------------------------------------
% Orientation-independent spin Hamiltonian components.
if ~higherZeemanPresent
  [H0,mux,muy,muz] = ham(Sys);
end

% Set up test field orientations
%-------------------------------------------------------------------------------
% Set of field vectors for the symmetry test operations. These directions are
% slightly off potential symmetry axes to detect rotation and reflection symmetries.
% Each column in FieldVecs tests a specific symmetry operation.
dtheta = 7.3234;  % small theta increment (degrees)
dphi = 6.3456;   % small phi increment (degrees)
theta = pi/180* [dtheta, dtheta, dtheta,dtheta/2,dtheta/2, dtheta,180-dtheta,180-dtheta,180-dtheta,90-dtheta, dtheta, dtheta];
phi = pi/180*([0, 120, 90, 90, 0, -90, 0, -90, 180, 90, -2*dphi, 180]+dphi);
Field = 350;  % field magnitude, mT
FieldVecs_S = Field*ang2vec(phi,theta);  % test field vectors, in (candidate) symmetry frame


% Determination of point group
%-------------------------------------------------------------------------------
% Start with C1 symmetry in molecular frame, then examine each one
% of the candidate symmetry frames for higher symmetry, and update.

highestSymmetry = 0;  % no symmetry (C1)
RMatrix = eye(3);

pointGroupNames = {'Ci','C2h','D2h','C4h','D4h','S6','D3d','C6h','D6h','Th','Oh','Dinfh','O3'};

for iFrame = 1:size(allFrames,3)  % loop over all candidate symmetry frame

  R_S2M = allFrames(:,:,iFrame);  % transformation matrix from symmetry to molecular frame
  
  % Transform field vectors from (candidate) symmetry frame to molecular frame
  B_M = R_S2M*FieldVecs_S;
  
  % Test for symmetry elements (compute eigenvalues lazily to avoid
  % expensive calls unless a particular branch requires them)
  eigsA = computeEigs(B_M(:,1));
  eigsB = computeEigs(B_M(:,2));
  eigsC = computeEigs(B_M(:,3));
  C3z = eqeig(eigsA,eigsB);  % C3 axis along z
  C4z = eqeig(eigsA,eigsC);  % C4 axis along z

  if ~C3z && ~C4z
    % Neither C3z nor C4z: Ci, C2h, D2h
    C2z = eqeig(eigsA, computeEigs(B_M(:,12)));  % C2 axis along z
    if ~C2z
      pointGroupId = 1;  % Ci
    else
      sigmaxz = eqeig(eigsA, computeEigs(B_M(:,11)));  % Mirror plane in xz plane
      if sigmaxz
        pointGroupId = 3;  % D2h
      else
        pointGroupId = 2;  % C2h
      end
    end

  elseif C3z && ~C4z
    % C3z, but no C4z: S6, D3d, Th, C6h, D6h
    sigmaxy = eqeig(eigsA, computeEigs(B_M(:,7)));  % Mirror plane in xy plane
    if sigmaxy
      % Th, C6h, D6h
      C2z = eqeig(eigsA, computeEigs(B_M(:,12)));  % C2 axis along z
      if C2z
        C2x = eqeig(eigsC, computeEigs(B_M(:,8)));  % C2 axis along x
        if C2x
          pointGroupId = 9;  % D6h
        else
          pointGroupId = 8;  % C6h
        end
      else
        pointGroupId = 10;   % Th
      end
    else
      % S6, D3d
      C2x = eqeig(eigsC, computeEigs(B_M(:,8)));  % C2 axis along x
      C2y = eqeig(eigsA, computeEigs(B_M(:,9)));  % C2 axis along y
      if C2x || C2y
        pointGroupId = 7;   % D3d
      else
        pointGroupId = 6;   % S6
      end
    end

  elseif ~C3z && C4z
    % C4z, but no C3z: C4h, D4h, Oh
    C2x = eqeig(eigsC, computeEigs(B_M(:,8)));  % C2 axis along x
    if C2x
      C3d = eqeig(computeEigs(B_M(:,10)), computeEigs(B_M([2 1 3],10)));  % C3 axis along [1;1;1]
      if C3d
        pointGroupId = 11;  % Oh
      else
        pointGroupId = 5;   % D4h
      end
    else
      pointGroupId = 4;     % C4h
    end

  else
    % Both C3z and C4z axes: Dinfh, O3
    Cinfx = eqeig(eigsC, computeEigs(B_M(:,4)));  % Cinf axis along x
    if Cinfx
      pointGroupId = 13;    % O3
    else
      pointGroupId = 12;    % Dinfh
    end
  end

  % Update if symmetry is higher than current best
  if pointGroupId>highestSymmetry
    highestSymmetry = pointGroupId;
    RMatrix = R_S2M;
  end
  
end


% Output assignment
%-------------------------------------------------------------------------------
% Assertion that a point group has been found.
if highestSymmetry==0
  error('No point group found! Save input spin system and report bug!');
end

Group = pointGroupNames{highestSymmetry};

  % Nested function: compute eigenvalues for a given field direction
  function eigvals = computeEigs(B)
    if higherZeemanPresent
      eigvals = eig(ham(Sys,B));
    else
      eigvals = eig(H0 - B(1)*mux - B(2)*muy - B(3)*muz);
    end
  end

end
%===============================================================================


%===============================================================================
% Determine if two sets of eigenvalues are (approximately) equal
function tf = eqeig(eigsA,eigsB)
threshold = 1e-12;
eigsA = sort(eigsA);
eigsB = sort(eigsB);
tf = norm(eigsA-eigsB) < threshold*(norm(eigsA)+norm(eigsB));
end


%===============================================================================
% Determine point group symmetry using geometric reasoning
function [symGrp,symFrame] = hamsymm_geom(Sys)

Sys = validatespinsys(Sys);
nElectrons = Sys.nElectrons;
nNuclei = Sys.nNuclei;
nElCouplings = nElectrons*(nElectrons-1)/2;
nNucCouplings = nNuclei*(nNuclei-1)/2;

% Determine symmetries of all tensors in spin system
%-------------------------------------------------------------------------------
% Preallocations
nTensors = nElectrons + sum(Sys.S>1/2) + nElCouplings + nNuclei + ...
  nNuclei*nElectrons + sum(Sys.I>1/2) + nNucCouplings;
tensorSym = zeros(1,nTensors);
tensorFrame = cell(1,nTensors);
tensorName = cell(1,nTensors);
idx = 1;

% g tensors
for iE = 1:nElectrons
  [tensorSym(idx),tensorFrame{idx}] = tensorsymmetry(Sys.g(iE,:),Sys.gFrame(iE,:));
  tensorName{idx} = sprintf('g%d',iE);
  idx = idx + 1;
end

% D tensors
for iE = 1:nElectrons
  if Sys.S(iE)>1/2
    [tensorSym(idx),tensorFrame{idx}] = tensorsymmetry(Sys.D(iE,:),Sys.DFrame(iE,:));
    tensorName{idx} = sprintf('D%d',iE);
    idx = idx + 1;
  end
end

% Electron-electron tensors
for iC = 1:nElCouplings
  [tensorSym(idx),tensorFrame{idx}] = tensorsymmetry(Sys.ee(iC,:),Sys.eeFrame(iC,:));
  tensorName{idx} = sprintf('ee%d',iC);
  idx = idx + 1;
end

% CSA tensors
for iN = 1:nNuclei
  [tensorSym(idx),tensorFrame{idx}] = tensorsymmetry(Sys.sigma(iN,:),Sys.sigmaFrame(iN,:));
  tensorName{idx} = sprintf('sigma%d',iN);
  idx = idx + 1;
end

% A tensors
for iN = 1:nNuclei
  eidx = 1:3;
  for iE = 1:nElectrons
    [tensorSym(idx),tensorFrame{idx}] = tensorsymmetry(Sys.A(iN,eidx),Sys.AFrame(iN,eidx));
    tensorName{idx} = sprintf('A%d%d',iE,iN);
    eidx = eidx + 3;
    idx = idx + 1;
  end
end

% Q tensors
for iN = 1:nNuclei
  if Sys.I(iN)>1/2
    [tensorSym(idx),tensorFrame{idx}] = tensorsymmetry(Sys.Q(iN,:),Sys.QFrame(iN,:));
    tensorName{idx} = sprintf('Q%d',iN);
    idx = idx + 1;
  end
end

% Nucleus-nucleus tensors
for iC = 1:nNucCouplings
  [tensorSym(idx),tensorFrame{idx}] = tensorsymmetry(Sys.nn(iC,:),Sys.nnFrame(iC,:));
  tensorName{idx} = sprintf('nn%d',iC);
  idx = idx + 1;
end

% Remove isotropic tensors
isotropic = tensorSym==0;
tensorSym(isotropic) = [];
tensorFrame(isotropic) = [];
tensorName(isotropic) = [];


% Compute total symmetry by iteratively pairwise combination of tensor symmetries
%-------------------------------------------------------------------------------
pointGroups = {'O3','Dinfh','D2h','C2h','Ci'};
grp = 0;
symFrame = eye(3);
for iTens = 1:numel(tensorSym)
  [grp,symFrame] = combinesymms(grp,symFrame,tensorSym(iTens),tensorFrame{iTens});
end
symGrp = pointGroups{grp+1};

end


%===============================================================================
% tensorsymmetry   Determine second-rank tensor symmetry
%
% [symGroup,symFrame] = tensorsymmetry(principalValues,EulerAngles)
%
% Input:
%   principalValue  3-element array of principal values
%   EulerAngles     3-element array of Euler angles, in radians
% Output:
%   symGroup        0 for isotropic/O3, 1 for axial/D2inf, 2 for rhombic/D2h)
%   symFrame        transformation matrix from symmetry frame to molecular frame

function [symGroup,symFrame] = tensorsymmetry(principalValues,EulerAngles)

O3 = 0;  % isotropic
Dinfh = 1;  % axial
D2h = 2;  % rhombic

if principalValues(1)==principalValues(2)
  if principalValues(1)==principalValues(3)
    symGroup = O3;
    xyz = [1 2 3];
  else
    symGroup = Dinfh;
    xyz = [1 2 3];
  end
elseif principalValues(1)==principalValues(3)
  symGroup = Dinfh;
  xyz = [2 3 1];
elseif principalValues(2)==principalValues(3)
  symGroup = Dinfh;
  xyz = [3 1 2];
else
  symGroup = D2h;
  xyz = [1 2 3];
end

R_M2S = erot(EulerAngles);
R_S2M = R_M2S.';
symFrame = R_S2M(:,xyz);

end


%===============================================================================
% combinesymms  Combine symmetries of two tensors
%
%      [combinedSym,combinedFrame] = combinesymms(Sym1,frame1,Sym2,frame2)
%
%      Computes the combined symmetry group of Sym1 and Sym2.
%      Sym1 and Sym2 are symmetry group IDs: 0 O3, 1 Dinfh, 2 D2h,
%      3 C2h, 4 Ci. frame1 and frame2 are the rotation matrices containing
%      the orientations of the principal symmetry axes in a reference frame
%      combinedSym is the combined symmetry group, and combinedFrame contains
%      the combined axes orientation. Sym2 must be either O3, Dinf or D2h.

function [combinedSym,combinedFrame] = combinesymms(Sym1,frame1,Sym2,frame2)

% Rotation matrix R = erot([alpha beta gamma])
%   Cols: tensor frame axes in reference frame representation
%   Rows: reference frame axes in tensor frame representation

% Symmetry group code abbreviations
O3 = 0; Dinfh = 1; D2h = 2; C2h = 3; Ci = 4;

% Angle thresholds for parallel and perpendicular tests
delta = 0.2;  % deg
parallelLimit = delta*pi/180;  % maximum
perpLimit = (90-delta)*pi/180;  % minimum

% Sort Sym1 and Sym2 so that Sym1 > Sym2
% (Sym1 is of lower symmetry than Sym2)
if Sym1<Sym2
  [Sym2,Sym1] = deal(Sym1,Sym2);
  [frame2,frame1] = deal(frame1,frame2);
end

% any + O3, Ci + any, O3 + any
%-------------------------------------------------------------------------------
% - If Sym1 is Ci, the total symmetry remains Ci, since it cannot be lower.
% - If Sym2 is O3, the total symmetry is unchanged.
if Sym1==Ci || Sym2==O3
  combinedSym = Sym1;
  combinedFrame = frame1;
  return
end
% - If Sym1 is O3, the total symmetry is the new one (either Dinfh or D2h).
if Sym1==O3
  combinedSym = Sym2;
  combinedFrame = frame2;
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
  
  z1 = frame1(:,3);
  z2 = frame2(:,3);
  z1z2Angle = real(acos(abs(z1'*z2)));
  if z1z2Angle<parallelLimit  % z axes parallel
    combinedSym = Dinfh;
    combinedFrame = frame1;
  elseif z1z2Angle>perpLimit  % z axes perpendicular
    combinedSym = D2h;
    newz = cross(z1,z2);
    newz = newz/norm(newz);
    combinedFrame = [cross(z1,newz) z1 newz];
  else % all other angles between z1 and z2
    combinedSym = C2h;
    newz = cross(z1,z2);
    newz = newz/norm(newz);
    combinedFrame = [cross(z1,newz) z1 newz];
  end

% [2] D2h + Dinfh
%-------------------------------------------------------------------------------
% - Dinfh z parallel to any D2h axis -> D2h as before
% - Dinfh z perpendicular to one D2h axis (i.e. Dinfh z in a D2h
%   mirror plane) -> C2h with the perpendicular D2h axis as z axis
% - otherwise Ci with arbitrary orientation
elseif Sym1==D2h && Sym2==Dinfh
  z2 = frame2(:,3); % Dinfh z axis
  z2Angles = real(acos(abs(z2'*frame1)));
  if any(z2Angles<parallelLimit)
    combinedSym = D2h;
    combinedFrame = frame1;
  elseif any(z2Angles>perpLimit) % z in a D2h sigma plane
    combinedSym = C2h;
    [~,newzidx] = max(z2Angles);
    idx = [2 3 1 2 3];
    combinedFrame = frame1(:,idx(newzidx+(0:2)));
  else
    combinedSym = Ci;
    combinedFrame = frame1; % can be arbitrary
  end
  
% [3] C2h + Dinfh
%-------------------------------------------------------------------------------
% - Dinfh z parallel to C2h z axis -> C2h as before
% - Dinfh z in the C2h mirror plane -> C2h as before
% - Dinfh z arbitrary -> Ci, arbitrary frame
elseif Sym1==C2h && Sym2==Dinfh
  
  z1 = frame1(:,3);
  z2 = frame2(:,3);
  z1z2Angle = real(acos(abs(z1'*z2)));
  if z1z2Angle<parallelLimit
    combinedSym = C2h;
    combinedFrame = frame1;
  elseif z1z2Angle>perpLimit
    combinedSym = C2h;
    combinedFrame = frame1;
  else
    combinedSym = Ci;
    combinedFrame = frame1; % can be arbitrary
  end
  
% [4] D2h + D2h
%-------------------------------------------------------------------------------
% - all axes collinear -> D2h as before
% - only one axis collinear -> C2h with collinar axis as z
% - all axes skew -> Ci, arbitrary frame
elseif Sym1==D2h && Sym2==D2h
  
  allAngles = real(acos(abs(frame2'*frame1)));
  % row 1: angles of x2, row 2: angles of y2, row 3: angles of z2 
  % col 1: angles of x1, col 2: angles of y1, col 3: angles of z1
  axisCollinear = min(allAngles)<parallelLimit;
  if all(axisCollinear)
    combinedSym = D2h;
    combinedFrame = frame1;
  elseif any(axisCollinear)
    combinedSym = C2h;
    % find collinear axes
    %[minAngles,minidx2] = min(allAngles);
    [minAngles] = min(allAngles);
    [~,newzidx1] = min(minAngles);
    %newzidx2 = minidx2(newzidx1);
    % Axes newzidx2 of frame 2 and newzidx1 of frame 1 are collinear
    idx = [2 3 1 2 3];
    combinedFrame = frame1(:,idx(newzidx1+(0:2)));
  else
    combinedSym = Ci;
    combinedFrame = frame1;% can be arbitrary
  end
  
% [5] C2h + D2h
%-------------------------------------------------------------------------------
% - any D2h axis collinar with C2h z axis -> C2h as before
% - otherwise Ci, arbitrary frame
elseif Sym1==C2h && Sym2==D2h

  z1 = frame1(:,3);
  z1Angles = real(acos(abs(z1'*frame2)));
  if any(z1Angles<parallelLimit)
    combinedSym = C2h;
    combinedFrame = frame1;
  else
    combinedSym = Ci;
    combinedFrame = frame1; % can be arbitrary
  end

else

  error('Cannot compute symmetry combination: %d %d',Sym1,Sym2);

end % Combination case switchyard

end

%===============================================================================
% Determine if the spin system contains isotropic terms only
function iso = isisotropic(Sys)

isiso = @(T) isfield(Sys,T) && any(any(diff(Sys.(T),[],2)));

iso = isiso('g') && isiso('D') && isiso('ee') && isiso('Q') && isiso('sigma') && isiso('nn');

if iso && isfield(Sys,'A')
  eidx = 1:3;
  for e = 1:Sys.nElectrons
    iso = iso && any(diff(Sys.A(:,eidx),[],2));
    eidx = eidx + 3;
  end
end

end

%===============================================================================
% Check any field stating with specific string have non-zero elements
function present = checkFields(Sys,fieldnamestart)

sysFields = fieldnames(Sys);
matchingFields = strncmp(sysFields,fieldnamestart,length(fieldnamestart)).';

present = false;
if ~any(matchingFields), return; end

for f = find(matchingFields)
  if any(Sys.(sysFields{f})(:))
    present = true;
    return
  end
end

end
