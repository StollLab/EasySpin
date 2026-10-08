% strains_setup  Set up strain data for a spin system
%
%   [StrainData,err] = strains_setup(SysIn,Sys)
%
%   Parses the strain fields (StrainPars, StrainFWHM, StrainCorr, StrainModes)
%   of the raw spin system SysIn and computes, for each strain parameter p_i,
%   the analytical derivative dT/dp_i of the affected interaction tensor T in
%   the molecular frame. Sys is the validated spin system (output of
%   validatespinsys).
%
%   Output:
%     StrainData.Q      n x m mode matrix (FWHM units), covariance = Q*Q.'
%     StrainData.Deriv  n-element struct array with fields
%                         type  'g','A','D','Q','ee','nn','sigma'
%                         idx   electron (g, D), [electron nucleus] (A),
%                               nucleus (Q, sigma), [e1 e2] (ee), [n1 n2] (nn)
%                         dT    3x3 derivative of the tensor in the molecular frame
%     err               error message, empty if no error

function [StrainData,err] = strains_setup(SysIn,Sys)

StrainData = struct('Q',[],'Deriv',[]);
err = '';

strainFields = {'StrainFWHM','StrainCorr','StrainModes'};
if ~isfield(SysIn,'StrainPars') || isempty(SysIn.StrainPars)
  for f = 1:numel(strainFields)
    if isfield(SysIn,strainFields{f})
      err = sprintf('Sys.%s is given, but Sys.StrainPars is missing or empty.',strainFields{f});
      return
    end
  end
  return
end

try
  [StrainData.Q,StrainData.Deriv] = setup(SysIn,Sys);
catch ME
  err = ME.message;
end

end


%-------------------------------------------------------------------------------
function [Q,Deriv] = setup(SysIn,Sys)

nElectrons = Sys.nElectrons;

% Nuclei as given by user (including spin-0 nuclei, which validatespinsys removes)
if isfield(SysIn,'Nucs') && ~isempty(SysIn.Nucs)
  rawNucList = nucstring2list(SysIn.Nucs);
else
  rawNucList = {};
end
rawI = nucdata(rawNucList);
nRawNuclei = numel(rawI);
keep = rawI~=0;
newNucIdx = cumsum(keep).*keep; % raw nucleus index -> validated nucleus index (0: removed)

refs = strains_parse(SysIn.StrainPars,SysIn,nElectrons,nRawNuclei);
nPars = numel(refs);
Q = strains_modes(SysIn,nPars);

elPairs = [];
if nElectrons>1, elPairs = nchoosek(1:nElectrons,2); end
rawNucPairs = [];
if nRawNuclei>1, rawNucPairs = nchoosek(1:nRawNuclei,2); end
nucPairs = [];
if Sys.nNuclei>1, nucPairs = nchoosek(1:Sys.nNuclei,2); end

Deriv = struct('type',cell(1,nPars),'idx',[],'dT',[]);
for p = 1:nPars
  ref = refs(p);
  f = ref.FieldName;
  r = ref.Subscripts(1);
  c = ref.Subscripts(2);
  isFrame = strcmp(ref.Form,'frame');
  if isFrame
    tensorName = f(1:end-5);
    instanceField = tensorName;
  else
    instanceField = f;
  end

  % Determine instance (electron, nucleus, or pair) and element
  switch instanceField
    case {'g','D'}
      type = instanceField;
      if any(strcmp(ref.Form,{'iso','D'})), e = ref.Index; k = 1; else, e = r; k = c; end
      idx = e;
      if strcmp(type,'D') && Sys.S(e)<1
        error('Sys.StrainPars: ''%s'' has no effect, since electron spin %d has S = 1/2.',ref.Name,e);
      end
    case {'Q','sigma'}
      type = instanceField;
      if any(strcmp(ref.Form,{'iso','eeqQ'})), n = ref.Index; k = 1; else, n = r; k = c; end
      n = mapnucleus(n,ref.Name);
      idx = n;
      if strcmp(type,'Q') && Sys.I(n)<1
        error('Sys.StrainPars: ''%s'' has no effect, since nucleus %s has I < 1.',ref.Name,Sys.Nucs{n});
      end
    case 'A'
      type = 'A';
      if isFrame
        n = r; e = ceil(c/3); k = c-3*(e-1);
      else
        switch ref.Form
          case 'iso1', n = c; e = 1; k = 1;
          case 'iso', n = r; e = c; k = 1;
          case 'axial', n = r; e = ceil(c/2); k = c-2*(e-1);
          case 'pv', n = r; e = ceil(c/3); k = c-3*(e-1);
          case 'sym', n = r; e = ceil(c/6); k = c-6*(e-1);
        end
      end
      n = mapnucleus(n,ref.Name);
      idx = [e n];
    case {'ee','J','dip','dvec'}
      type = 'ee';
      if any(strcmp(ref.Form,{'iso','J','dip1'})), pair = ref.Index; k = 1; else, pair = r; k = c; end
      idx = elPairs(pair,:);
    case 'nn'
      type = 'nn';
      if strcmp(ref.Form,'iso'), rawPair = ref.Index; k = 1; else, rawPair = r; k = c; end
      n1 = mapnucleus(rawNucPairs(rawPair,1),ref.Name);
      n2 = mapnucleus(rawNucPairs(rawPair,2),ref.Name);
      pair = find(nucPairs(:,1)==n1 & nucPairs(:,2)==n2);
      idx = [n1 n2];
  end

  % Euler angles of the tensor frame (zero for full tensors)
  switch type
    case {'g','D'}, angles = Sys.([type 'Frame'])(e,:);
    case {'Q','sigma'}, angles = Sys.([type 'Frame'])(n,:);
    case 'A', idxE = 3*(e-1)+(1:3); angles = Sys.AFrame(n,idxE);
    case {'ee','nn'}, angles = Sys.([type 'Frame'])(pair,:);
  end

  if isFrame
    % Tensor in its own frame
    switch type
      case {'g','D'}, P = diag(Sys.(type)(e,:));
      case {'Q','sigma'}, P = diag(Sys.(type)(n,:));
      case 'A', P = diag(Sys.A(n,idxE));
      case 'nn', P = diag(Sys.nn(pair,:));
      case 'ee'
        if Sys.fullee
          P = Sys.ee(3*(pair-1)+(1:3),:);
        else
          P = diag(Sys.ee(pair,:));
        end
    end
    % Derivative with respect to Euler angle k: T = R*P*R.'
    [R,Rdot] = rotderiv(angles,k);
    dT = Rdot*P*R.' + R*P*Rdot.';
    scale = norm(P,'fro');
  else
    if strcmp(type,'Q'), I = Sys.I(n); else, I = []; end
    R = erot(angles).'; % tensor frame -> molecular frame
    dP = paramderiv(ref,k,SysIn,I);
    dT = R*dP*R.';
    scale = 1;
  end

  if norm(dT,'fro') <= 1e-10*max(scale,eps)
    error('Sys.StrainPars: ''%s'' has no effect on the spin Hamiltonian.',ref.Name);
  end

  Deriv(p).type = type;
  Deriv(p).idx = idx;
  Deriv(p).dT = dT;
end

  %-----------------------------------------------------------------------------
  function n = mapnucleus(nRaw,name)
  n = newNucIdx(nRaw);
  if n==0
    error('Sys.StrainPars: ''%s'' refers to a nucleus with spin 0.',name);
  end
  if Sys.n(n)>1
    error('Sys.StrainPars: ''%s'' refers to a set of equivalent nuclei (Sys.n>1), which is not supported.',name);
  end
  end

end


%-------------------------------------------------------------------------------
% Derivative of the tensor in its own frame with respect to a tensor value
function dP = paramderiv(ref,k,SysIn,I)

Pi = @(j) full(sparse(j,j,1,3,3));
switch ref.Form
  case {'iso','iso1','J'}
    dP = eye(3);
  case 'axial'
    if k==1, dP = diag([1 1 0]); else, dP = diag([0 0 1]); end
  case 'pv'
    dP = Pi(k);
  case 'sym'
    pairs = [1 1; 2 2; 3 3; 1 2; 1 3; 2 3];
    dP = zeros(3);
    dP(pairs(k,1),pairs(k,2)) = 1;
    dP(pairs(k,2),pairs(k,1)) = 1;
  case 'D'
    dP = diag([-1 -1 2]/3);
  case 'DE'
    if k==1, dP = diag([-1 -1 2]/3); else, dP = diag([1 -1 0]); end
  case {'eeqQ','eeqQeta'}
    % Q principal values: eeqQ/(4I(2I-1))*[-1+eta, -1-eta, 2]
    if strcmp(ref.Form,'eeqQ')
      eeqQ = SysIn.Q(ref.Index);
      eta = 0;
    else
      n = ref.Subscripts(1);
      eeqQ = SysIn.Q(n,1);
      eta = SysIn.Q(n,2);
    end
    pre = 1/(4*I*(2*I-1));
    if k==1
      dP = pre*diag([-1+eta, -1-eta, 2]);
    else
      dP = pre*eeqQ*diag([1 -1 0]);
    end
  case 'dip1'
    dP = diag([1 1 -2]);
  case 'dip2'
    if k==1, dP = diag([1 1 -2]); else, dP = diag([1 -1 0]); end
  case 'dip3'
    dP = Pi(k) - eye(3)/3;
  case 'dvec'
    % ee = J*eye(3) + [0 d3 -d2; -d3 0 d1; d2 -d1 0] + diag(dip)
    switch k
      case 1, dP = [0 0 0; 0 0 1; 0 -1 0];
      case 2, dP = [0 0 -1; 0 0 0; 1 0 0];
      case 3, dP = [0 1 0; -1 0 0; 0 0 0];
    end
  otherwise
    error('Sys.StrainPars: Sys.%s has an unsupported size.',ref.FieldName);
end

end


%-------------------------------------------------------------------------------
% Rotation matrix R = erot(angles).' (tensor frame -> molecular frame) and its
% derivative with respect to Euler angle k (1: alpha, 2: beta, 3: gamma).
% erot = Rg*Rb*Ra, see erot.m
function [R,Rdot] = rotderiv(angles,k)

ca = cos(angles(1)); sa = sin(angles(1));
cb = cos(angles(2)); sb = sin(angles(2));
cg = cos(angles(3)); sg = sin(angles(3));

Ra = [ca sa 0; -sa ca 0; 0 0 1];
Rb = [cb 0 -sb; 0 1 0; sb 0 cb];
Rg = [cg sg 0; -sg cg 0; 0 0 1];

switch k
  case 1, dRa = [-sa ca 0; -ca -sa 0; 0 0 0]; dE = Rg*Rb*dRa;
  case 2, dRb = [-sb 0 -cb; 0 0 0; cb 0 -sb]; dE = Rg*dRb*Ra;
  case 3, dRg = [-sg cg 0; -cg -sg 0; 0 0 0]; dE = dRg*Rb*Ra;
end

R = (Rg*Rb*Ra).';
Rdot = dE.';

end
