% validatespinsys   Validation of spin system structure
%
%   [FullSys,err] = validatespinsys(Sys)
%
%   Returns a non-empty error string in err if spin system Sys is not valid.
%   FullSys is the processed spin system. All missing optional fields
%   are supplemented, and several other fields are added, such as
%
%     nElectrons, nNuclei, Spins, nStates
%     fullg, fullA, fullD, fullQ, fullnn
%     I, gn

function [FullSys,err] = validatespinsys(Sys)

FullSys = [];
err = '';

if ~isstruct(Sys) || (numel(Sys)~=1)
  err = 'Spin system must be a structure.';
  return
end

% whether Sys is being reprocessed after removal of nuclei
reprocessing = false;
if isfield(Sys,'processed')
  if Sys.processed
    FullSys = Sys;
    return
  else
    reprocessing = true;
  end
end

% Mex compilation check
%-------------------------------------------------------------------------------
fileName = 'cubicsolve';
if exist(fileName,'file')~=3
 easyspincompile;
  if exist(fileName,'file')~=3
    error('EasySpin: Generation of mex files failed.');
  end
end

% Spell check field names (capitalization)
%-------------------------------------------------------------------------------
correctFields = {'S','Nucs','Abund','n',...
  'g','g_','gFrame','gStrain',...
  'D','DFrame','DStrain',...
  'ee','J','dip','dvec','ee2','eeFrame',...
  'A','A_','AFrame','AStrain',...
  'Q','QFrame',...
  'HStrain',...
  'L', 'soc', 'orf',...
  'lw','lwpp','lwEndor',...
  'tcorr','logtcorr','Diff','logDiff'};
fieldlist = @(str,irange)arrayfun(@(x)sprintf('%s%d',str,x),irange,'UniformOutput',false);
correctFields = [correctFields fieldlist('B',1:12)];
correctFields = [correctFields fieldlist('CF',1:12)];

givenFields = fieldnames(Sys);
for f = 1:numel(givenFields)
  givField = givenFields{f};
  if strcmp(givField,'ZeemanFreq')
    err = 'Field Sys.ZeemanFreq can only be used in conjunction with the function spidyan.';
    return
  end
  idx = find(strcmpi(givField,correctFields));
  % check if there is a case-insensitive match
  if idx
    corrField = correctFields{idx};
    if ~strcmp(givField,corrField)
      % Wrong capitalization
      error('Fix capitalization: Sys.%s should be Sys.%s',givField,corrField);
    else
      % Correct capitalization
    end
  end
end

for ind = find((strncmpi(givenFields,'Ham',3)))
  if isempty(ind), break; end
  field = givenFields{ind};
  if length(field)~= 6 
    if str2double(field(4))+str2double(field(5))<10  
      err = sprintf('Wrong length of Sys.%s entry, should be Hamxyz (with x,y,z integer numbers)',field);
      return
    else
      if length(field)~= 7
        err = sprintf('Wrong length of Sys.%s entry, should be Hamxyz (with x,y,z integer numbers)',field);
        return
      end
    end
  end
  if ~strncmp(field,'Ham',3)
      % Wrong capitalization
      err = sprintf('Fix capitalization: Sys.%s should be Sys.%s',field,['Ham', field(4:end)]);
      return
  end
end


% Electron spins
%-------------------------------------------------------------------------------
% If S is missing, set it to 1/2
if ~isfield(Sys,'S')
  Sys.S = 1/2;
end

% Guard against invalid type
if isempty(Sys.S) || any(~isreal(Sys.S)) || any(mod(real(Sys.S),1/2)) || any(Sys.S<0)
  err = 'Electron spins in S must be positive integer multiples of 1/2.';
  return
end
if any(Sys.S==0)
  err = 'Systems with zero spin are not allowed.';
  return
end

nElectrons = numel(Sys.S);
Sys.nElectrons = nElectrons;


% g tensor(s) (Sys.g, Sys.g_, Sys.gFrame)
%-------------------------------------------------------------------------------

if isfield(Sys,'g')
  
  if isfield(Sys,'g_')
    err = 'Sys.g and Sys.g_ are given. Remove one of them.';
    return
  end
  
  Sys.fullg = issize(Sys.g,[3*nElectrons 3]);
  if Sys.fullg
    % full g tensors
  elseif numel(Sys.g)==nElectrons
    % isotropic g factors
    Sys.g = Sys.g(:)*[1 1 1];
  elseif issize(Sys.g,[nElectrons 2])
    % axial tensors
    Sys.g = Sys.g(:,[1 1 2]);
  elseif issize(Sys.g,[nElectrons 3])
    % orthorhombic tensors
  else
    err = 'Sys.g has wrong size.';
    return
  end

elseif isfield(Sys,'g_')

  Sys.fullg = false;
  if numel(Sys.g_)==nElectrons
    % isotropic g factors
    Sys.g_(nElectrons,3) = 0;
  elseif issize(Sys.g_,[nElectrons 2])
    % axial tensors
    Sys.g_(nElectrons,3) = 0;
  elseif issize(Sys.g_,[nElectrons 3])
    % orthorhombic tensors
  else
    err = ('Sys.g_ has wrong size.');
    return
  end
  for iElectron = 1:nElectrons
    g_spherical = Sys.g_(iElectron,:);
    g_cartesian = g_spherical(1) + ...
                  g_spherical(2)*[-1 -1 2] + ...
                  g_spherical(3)*[-1 +1 0];
    Sys.g(iElectron,:) = g_cartesian;
  end
  
else
  
  % Supplement g
  Sys.fullg = false;
  if any(strncmp(fieldnames(Sys),'Ham',3))
    Sys.g = zeros(nElectrons,3);
  else
    Sys.g = gfree*ones(nElectrons,3);
  end
  
end

% Euler angles for g tensor(s)
err = pa_obsolete_message(Sys,'gpa','gFrame');
if ~isempty(err); return; end
if ~isfield(Sys,'gFrame') || isempty(Sys.gFrame)
  Sys.gFrame = zeros(nElectrons,3);
end
err = sizecheck(Sys,'gFrame',[nElectrons 3]);
if ~isempty(err); return; end


% Zero-field splittings (Sys.D, Sys.DFrame)
%-------------------------------------------------------------------------------

% Supplement partial or missing D
if ~isfield(Sys,'D')
  Sys.D = zeros(nElectrons,3);
  Sys.fullD = false;
else
  Sys.fullD = issize(Sys.D,[3*nElectrons 3]);
  if Sys.fullD
    % full D tensors
  elseif numel(Sys.D)==nElectrons
    % only D values
    Sys.D = Sys.D(:)*[-1/3,-1/3,+2/3];
  elseif issize(Sys.D,[nElectrons 2])
    % D and E values
    Sys.D = Sys.D(:,1)*[-1/3,-1/3,+2/3] + Sys.D(:,2)*[+1,-1,0];
  elseif issize(Sys.D,[nElectrons 3])
    % D principal values
  else
    err = ('Sys.D has wrong size.');
    return
  end
end

% Euler angles for D tensor(s)
err = pa_obsolete_message(Sys,'Dpa','DFrame');
if ~isempty(err); return; end
if ~isfield(Sys,'DFrame') || isempty(Sys.DFrame)
  Sys.DFrame = zeros(nElectrons,3);
end
err = sizecheck(Sys,'DFrame',[nElectrons 3]);
if ~isempty(err); return; end


% High-order zero-field terms (Sys.B*)
%-------------------------------------------------------------------------------

if isfield(Sys,'aF') || isfield(Sys,'aFFrame')
  error('Sys.aF and Sys.aFFrame are no longer supported. Use Sys.B4 and Sys.B4Frame instead.');
end
if isfield(Sys,'BFrame') && ~reprocessing
  err = 'Sys.BFrame is not supported. Use Sys.B1Frame, Sys.B2Fame, etc instead.';
  return
end

% B1, B2, B3, etc.
Sys.B = [];
D_present = any(Sys.D(:));
for k = 1:12
  fieldname = sprintf('B%d',k);
  if ~isfield(Sys,fieldname), continue; end
  Bk = Sys.(fieldname);
  
  if k==2 && any(Bk(:)) && D_present
    err = 'Cannot use Sys.D and Sys.B2 simultaneously. Remove one of them.';
    return
  end
  
  if size(Bk,1)~=nElectrons
    err = sprintf('Field Sys.%s has to have %d rows, since there are %d electron spins.',fieldname,nElectrons,nElectrons);
    return
  end
  
  if size(Bk,2)==1
    Bk = [zeros(nElectrons,k) Bk(:) zeros(nElectrons,k)];
  elseif size(Bk,2)==2*k+1
    % full form
  else
    err = sprintf('Field Sys.%s has %d instead of %d columns.',fieldname,size(Bk,2),2*k+1);
    return
  end
  
  if any(~isreal(Bk))
    err = sprintf('Field Sys.%s contains complex numbers. Only real ones are possible',fieldname);
    return
  end
  
  Sys.(fieldname) = Bk;
  Sys.B{k} = Bk;
  
  fieldname = sprintf('B%dFrame',k);
  if isfield(Sys,fieldname) && any(Sys.(fieldname)~=0)
    Sys.BFrame{k} = Sys.(fieldname);
  else
    Sys.BFrame{k} = [0 0 0];
  end
    
end


% Electron-electron coouplings (Sys.ee, Sys.J, Sys.dvec, Sys.dip, Sys.eeFrame)
%-------------------------------------------------------------------------------
if ~isfield(Sys,'fullee'), Sys.fullee = false; end
if nElectrons>1 && ~reprocessing
  
  eeMatrix = isfield(Sys,'ee');
  JdD = (isfield(Sys,'J') && ~isempty(Sys.J)) || ...
        (isfield(Sys,'dvec') && ~isempty(Sys.dvec)) || ...
        (isfield(Sys,'dip') && ~isempty(Sys.dip));
  
  if ~eeMatrix && ~JdD
    err = 'Spin system contains 2 or more electron spins, but coupling terms are missing (ee; or J, dip, dvec)!';
    return
  end
  
  if eeMatrix && JdD
    err = 'Both Sys.ee and (Sys.J,Sys.dip,Sys.dvec) are given - use only one or the other!';
    return
  end
  
  nElPairs = nElectrons*(nElectrons-1)/2;
  
  if eeMatrix
    % Bilinear coupling defined via Sys.ee
    
    % Expand isotropic couplings into 3 equal principal values
    if numel(Sys.ee)==nElPairs
      Sys.ee = Sys.ee(:)*[1 1 1];
    end
    
    fullee = issize(Sys.ee,[3*nElPairs,3]);
    Sys.fullee = fullee;
    if ~fullee
      err = sizecheck(Sys,'ee',[nElPairs 3]);
      if ~isempty(err), return; end
    end
    
    err = pa_obsolete_message(Sys,'eepa','eeFrame');
    if ~isempty(err), return; end
    
  else
    % Bilinear coupling defined via J, dip, and dvec
    % J:    isotropic exchange +J*S1*S2
    % dip:  dipolar coupling
    %        - 1 value: axial component
    %        - 2 values: axial and rhombic component
    %        - 3 values: principal values of dipolar tensor
    % dvec: antisymmetric exchange dvec.(S1xS2)
    
    % Size check on list of isotropic exchange coupling constants
    if ~isfield(Sys,'J'), Sys.J = zeros(1,nElPairs); end
    Sys.J = Sys.J(:);
    err = sizecheck(Sys,'J',[nElPairs 1]);
    if ~isempty(err), return; end
    
    % Size check on list of antisymmetric exchange vectors
    if ~isfield(Sys,'dvec'), Sys.dvec = zeros(nElPairs,3); end
    err = sizecheck(Sys,'dvec',[nElPairs,3]);
    if ~isempty(err), return; end

    % Size check on dipolar tensor diagonals
    if ~isfield(Sys,'dip'), Sys.dip = zeros(nElPairs,3); end
    if numel(Sys.dip)==nElPairs
      Sys.dip = Sys.dip(:);
    end
    if size(Sys.dip,1)~=nElPairs
      err = sprintf('Sys.dip must contain %d rows, since there are %d unique pairs of electron spins.',nElPairs,nElPairs);
      return
    end

    if isfield(Sys,'eeD')
      err = 'Sys.eeD is obsolete. Use Sys.dip instead.';
      return
    end
    
    % Convert axial/rhombic components to principal values
    switch size(Sys.dip,2)
      case 1
        Sys.dip = Sys.dip*[1 1 -2];
      case 2
        Sys.dip = Sys.dip(:,1)*[1 1 -2] + Sys.dip(:,2)*[+1 -1 0];
      case 3
        % Remove isotropic component to guarantee zero traces of dipolar tensors
        Sys.dip = Sys.dip - repmat(mean(Sys.dip,2),1,3);
      otherwise
        err = 'Sys.dip must contain 1, 2, or 3 columns.';
        return
    end
    
    % Combine (Sys.J,Sys.dip,Sys.dvec) into full interaction matrix in Sys.ee
    fullee = any(Sys.dvec(:)~=0);
    Sys.fullee = fullee;
    if fullee
      idx = 1:3;
      for iPair = 1:nElPairs
        J = Sys.J(iPair);
        d = Sys.dvec(iPair,:);
        ee = J*eye(3) + ...
          [0 d(3) -d(2); -d(3) 0 d(1); d(2) -d(1) 0] + ...
          diag(Sys.dip(iPair,:));
        Sys.ee(idx,:) = ee;
        idx = idx + 3;
      end
    else
      for iPair = 1:nElPairs
        Sys.ee(iPair,:) = Sys.J(iPair) + Sys.dip(iPair,:);
      end
    end
    
  end
  
  % Check for eeFrame, and supplement or error if necessary
  if fullee
    if isfield(Sys,'eeFrame')
      err = sprintf('Full matrices are specified in ee, so eeFrame is not allowed.');
      if ~isempty(err), return; end
    end
  else
    if ~isfield(Sys,'eeFrame'), Sys.eeFrame = zeros(nElPairs,3); end
    err = sizecheck(Sys,'eeFrame',[nElPairs 3]);
    if ~isempty(err), return; end
  end
end


% Isotropic biquadratic exchange (Sys.ee2)
%-------------------------------------------------------------------------------
if nElectrons>1
  
  nElPairs = nElectrons*(nElectrons-1)/2;
  if ~isfield(Sys,'ee2')
    Sys.ee2 = zeros(nElPairs,1);
  else
    Sys.ee2 = Sys.ee2(:);
  end
  err = sizecheck(Sys,'ee2',[nElPairs 1]);
  if ~isempty(err), return; end
  
end


% Nuclear spins (Sys.Nucs, Sys.n, Sys.gnscale)
%===============================================================================

if ~isfield(Sys,'Nucs')
  Sys.Nucs = '';
end

if isempty(Sys.Nucs)
  if isfield(Sys,'A')   && ~isempty(Sys.A) || ...
     isfield(Sys,'AFrame') && ~isempty(Sys.AFrame) || ...
     isfield(Sys,'Q')   && ~isempty(Sys.q) || ...
     isfield(Sys,'QFrame') && ~isempty(Sys.QFrame)
    err = 'The system contains A and/or Q fields, but no nucleus is specified!';
    if ~isempty(err), return; end
  end
end

if ~isempty(Sys.Nucs)
  Sys.Nucs = nucstring2list(Sys.Nucs);
end

[I,gn] = nucdata(Sys.Nucs);

Sys.I = I;
Sys.gn = gn;
nNuclei = numel(I);
Sys.nNuclei = nNuclei;

if isfield(Sys,'n')
  if isempty(Sys.n), Sys.n = ones(1,nNuclei); end
  err = sizecheck(Sys,'n',[1 nNuclei]);
  if ~isempty(err), return; end
else
  Sys.n = ones(1,nNuclei);
end

if isfield(Sys,'gnscale')
  if numel(Sys.gnscale)<nNuclei
    err = ('Incorrect number of elements in Sys.gnscale.');
    if ~isempty(err), return; end
  end
else
  if nNuclei>0
    Sys.gnscale = ones(1,nNuclei);
  else
    Sys.gnscale = [];
  end
end


% Chemical shielding tensor (Sys.sigma, Sys.sigmaFrame)
%-------------------------------------------------------------------------------
if isfield(Sys,'sigma')
    
  Sys.fullsigma = issize(Sys.sigma,[3*nNuclei 3]);
  if Sys.fullsigma
    % full CS tensors
  elseif numel(Sys.sigma)==nNuclei
    % isotropic CS tensors
    Sys.sigma = Sys.sigma(:)*[1 1 1];
  elseif issize(Sys.sigma,[nNuclei 2])
    % axial tensors
    Sys.sigma = Sys.sigma(:,[1 1 2]);
  elseif issize(Sys.sigma,[nNuclei 3])
    % orthorhombic tensors
  else
    err = 'Sys.sigma has wrong size.';
    return
  end

else
  
  % Supplement Sys.sigma
  Sys.fullsigma = false;
  Sys.sigma = ones(nNuclei,3);
  
end

if ~isfield(Sys,'sigmaFrame') || isempty(Sys.sigmaFrame)
  Sys.sigmaFrame = zeros(nNuclei,3);
end
err = sizecheck(Sys,'sigmaFrame',[nNuclei 3]);
if ~isempty(err); return; end


% Hyperfine couplings (Sys.A, Sys.A_, Sys.AFrame)
%-------------------------------------------------------------------------------
Sys.fullA = false;
if nNuclei>0
  
  if ~isfield(Sys,'A') && ~isfield(Sys,'A_')
    err = sprintf('No hyperfine tensors A given for the %d electron spins and %d nuclei in the system!',nElectrons, nNuclei);
    if ~isempty(err), return; end
  end
  
  if isfield(Sys,'A') && isfield(Sys,'A_')
    err = sprintf('Hyperfine data given in both Sys.A and Sys.A_. Please remove one of them.');
    if ~isempty(err), return; end
  end

  % Spherical representation  [aiso T rho]
  if isfield(Sys,'A_')
    
    if issize(Sys.A_,[1 nNuclei])
      % Allow simple one-row syntax in the case of 1 eletron spin
      Sys.A_ = Sys.A_.';
      Sys.A_(:,3) = 0;
    elseif issize(Sys.A_,[nNuclei,nElectrons])
      % Expand aiso to [aiso 0 0]
      A_ = Sys.A_;
      idx = 1;
      for iElectron=1:nElectrons
        Sys.A_(:,idx) = A_(:,iElectron);
        Sys.A_(:,idx+2) = 0;
        idx = idx + 3;
      end
    elseif issize(Sys.A_,[nNuclei,2*nElectrons])
      % Expand [aiso T] into [aiso T 0]
      A_ = Sys.A_;
      idx1 = 1;
      idx2 = 1;
      for iElectron=1:nElectrons
        Sys.A_(:,idx1) = A_(:,idx2);
        Sys.A_(:,idx1+1) = A_(:,idx2+1);
        Sys.A_(:,idx1+2) = 0;
        idx1 = idx1 + 3;
        idx2 = idx2 + 2;
      end
    elseif issize(Sys.A_,[nNuclei,3*nElectrons])
      % [aiso T rho]
    else
      err = ('Wrong size of the A_ hyperfine array in the spin system.');
      if ~isempty(err), return; end
    end
    
    % Convert to cartesian
    idx = 1:3;
    for iElectron = 1:nElectrons
      for iNucleus = 1:nNuclei
        A_spherical = Sys.A_(iNucleus,idx);
        A_cartesian = A_spherical(1) + [-1 -1 2]*A_spherical(2) + [-1 +1 0]*A_spherical(3);
        Sys.A(iNucleus,idx) = A_cartesian;
      end
      idx = idx + 3;
    end
    
  else

    % Cartesian representation  [Ax Ay Az]
    if ~isnumeric(Sys.A)
      err = 'Sys.A must be a numeric array.';
      return
    end
    
    if issize(Sys.A,[3*nNuclei,3*nElectrons])
      % Full A matrices
      Sys.fullA = true;
    elseif issize(Sys.A,[1 nNuclei])
      % Allow simple one-row syntax in the case of 1 eletron spin
      if nElectrons==1
        Sys.A = Sys.A(:)*[1 1 1];
      else
        err = 'Size of Sys.A matrix is inconsistent with number of electrons and nuclei.';
        if ~isempty(err), return; end
      end
    elseif issize(Sys.A,[nNuclei,nElectrons])
      % Expand isotropic A into 3 equal principal values
      Sys.A = kron(Sys.A,[1 1 1]);
    elseif issize(Sys.A,[nNuclei,2*nElectrons])
      % Expand axial A into 3 principal values
      idx = [1 1 2];
      for k = 2:nElectrons
        idx = [idx, 2*k-[1 1 0]];
      end
      Sys.A = Sys.A(:,idx);
    elseif issize(Sys.A,[nNuclei,3*nElectrons])
      % Three principal values for each A tensor given
    else
      err = sprintf('Size of Sys.A (%dx%d) is inconsistent with number of nuclei (%d) and electrons (%d).',size(Sys.A,1),size(Sys.A,2),nNuclei,nElectrons);
      if ~isempty(err), return; end
    end
    
  end
  
  % Euler angles for A tensor(s)
  err = pa_obsolete_message(Sys,'Apa','AFrame');
  if ~isempty(err); return; end
  if ~isfield(Sys,'AFrame') || isempty(Sys.AFrame)
    Sys.AFrame = zeros(nNuclei,3*nElectrons);
  end
  err = sizecheck(Sys,'AFrame',[nNuclei,3*nElectrons]);
  if ~isempty(err), return; end
  
end


% Nuclear quadrupole interaction (Sys.Q, Sys.QFrame)
%-------------------------------------------------------------------------------
Sys.fullQ = false;
if nNuclei>0
  
  if ~isfield(Sys,'Q')
    
    Sys.Q = zeros(nNuclei,3);
    Sys.fullQ = false;
    
  else
    
    Sys.fullQ = issize(Sys.Q,[3*nNuclei 3]);
    if ~Sys.fullQ
      % Supplement eta=0 if not given
      if numel(Sys.Q)==nNuclei
        Sys.Q = Sys.Q(:);
        Sys.Q(:,2) = 0;
      end
      % Convert eeqQ,eta -> principal values
      if issize(Sys.Q,[nNuclei,2])
        for iN = 1:nNuclei
          I = Sys.I(iN);
          if I<1
            Sys.Q(iN,1:3) = 0;
          else
            eeqQh = Sys.Q(iN,1); % eeqQ/h, in MHz
            eta = Sys.Q(iN,2);
            Sys.Q(iN,1:3) = eeqQh/(4*I*(2*I-1)) * [-1+eta, -1-eta, 2];
          end
        end
      end
      
      err = sizecheck(Sys,'Q',[nNuclei,3]);
      if ~isempty(err), return; end
    end
    
  end
  
  % Assert Q matrix is symmetric
  if Sys.fullQ
    for iNuc = 1:nNuclei
      Q_ = Sys.Q(3*(iNuc-1)+(1:3),:);
      if norm(Q_-Q_.')/norm(Q_)>1e-5
        err = 'Sys.Q contains asymmetric full Q matrix. Only symmetric Q matrices are allowed.';
        if ~isempty(err), return; end
      end
    end
  end
  
  % Euler angles for Q tensor(s)
  err = pa_obsolete_message(Sys,'Qpa','QFrame');
  if ~isempty(err); return; end
  if ~isfield(Sys,'QFrame') || isempty(Sys.QFrame)
    Sys.QFrame = zeros(nNuclei,3);
  end
  err = sizecheck(Sys,'QFrame',[nNuclei 3]);
  if ~isempty(err); return; end
  
end


% Nuclear-nuclear couplings (Sys.nn, Sys.nnFrame)
%-------------------------------------------------------------------------------
Sys.fullnn = false;
if nNuclei<2
  
  if isfield(Sys,'nn') && ~isempty(Sys.nn) && any(Sys.nn(:))
    err = 'Nuclear-nuclear couplings specified in Sys.nn, but fewer than two nuclei given.';
    return
  end
  
else
  
  % Bilinear coupling defined via Sys.nn
  nNucPairs = nNuclei*(nNuclei-1)/2;
  
  if isfield(Sys,'nn') && ~isempty(Sys.nn) && any(Sys.nn(:))
    
    % Expand isotropic couplings into 3 equal principal values
    if numel(Sys.nn)==nNucPairs
      Sys.nn = Sys.nn(:)*[1 1 1];
    end
    
    % Size checks for Sys.nn
    Sys.fullnn = issize(Sys.nn,[3*nNucPairs,3]);
    if ~Sys.fullnn
      err = sizecheck(Sys,'nn',[nNucPairs 3]);
      if ~isempty(err), return; end
    end
    
  else
    Sys.fullnn = false;
    Sys.nn = zeros(nNucPairs,3);
    Sys.nnFrame = zeros(nNucPairs,3);
  end
  
  % Check for nnFrame, supplement or error if necessary
  if Sys.fullnn
    if isfield(Sys,'nnFrame') && ~isempty(Sys.nnFrame)
      err = sprintf('Full matrices are specified in Sys.nn, so nnFrame is not allowed.');
      if ~isempty(err), return; end
    end
  else
    if ~isfield(Sys,'nnFrame'), Sys.nnFrame = zeros(nNucPairs,3); end
    err = sizecheck(Sys,'nnFrame',[nNucPairs 3]);
    if ~isempty(err), return; end
  end
  
end


% Remove spin-zero nuclei
%-------------------------------------------------------------------------------
rmv = Sys.I==0;

if any(rmv)

  Sys.Nucs(rmv) = [];
  Sys.I(rmv) = [];
  Sys.gn(rmv) = [];
  Sys.nNuclei = numel(Sys.gn);
  Sys.gnscale(rmv) = [];
  Sys.n(rmv) = [];
  
  if Sys.fullA
    rmvfull = logical(kron(rmv(:),true(3,1)));
    Sys.A(rmvfull,:) = [];
  else
    Sys.A(rmv,:) = [];
    Sys.AFrame(rmv,:) = [];
  end
    
  if Sys.fullQ
    rmvfull = logical(kron(rmv(:),true(3,1)));
    Sys.Q(rmvfull,:) = [];
  else
    Sys.Q(rmv,:) = [];
    Sys.QFrame(rmv,:) = [];
  end
  
end


% Broadenings (strains and convolution line widths)
%===============================================================================

% Check Sys.lw
if ~isfield(Sys,'lw'), Sys.lw = [0 0]; end
if numel(Sys.lw)==1, Sys.lw(2) = 0; end
if numel(Sys.lw)~=2, err = ('System.lw has wrong size.'); end
if ~isempty(err), return; end
if any(Sys.lw<0), err = ('System.lw cannot be negative.'); end
if ~isempty(err), return; end

% Check Sys.lwEndor
if ~isfield(Sys,'lwEndor'), Sys.lwEndor = [0 0]; end
if numel(Sys.lwEndor)==1, Sys.lwEndor(2) = 0; end
if numel(Sys.lwEndor)~=2, err = ('System.lwEndor has wrong size.'); end
if ~isempty(err), return; end
if any(Sys.lwEndor<0), err = ('System.lwEndor cannot be negative.'); end
if ~isempty(err), return; end

% Check Sys.lwpp
if ~isfield(Sys,'lwpp'), Sys.lwpp = [0 0]; end
if numel(Sys.lwpp)==1, Sys.lwpp(2) = 0; end
if numel(Sys.lwpp)~=2, err = ('System.lwpp has wrong size.'); end
if ~isempty(err), return; end
if any(Sys.lwpp<0), err = ('Linewidths cannot be negative.'); end
if ~isempty(err), return; end

% Convert Sys.lwpp to Sys.lw (the latter is used internally)
if any(Sys.lwpp)
  if any(Sys.lw)
    err = ('System.lw and Sys.lwpp cannot be used simultaneously.');
    if ~isempty(err), return; end
  end
  Sys.lw = Sys.lwpp .* [sqrt(2*log(2)), sqrt(3)];
  Sys.lwpp = [0 0];
end


% g strain (Sys.gStrain)
%-------------------------------------------------------------------------------
if ~isfield(Sys,'gStrain') || isempty(Sys.gStrain)
  Sys.gStrain = zeros(nElectrons,3);
end

[n1,n2] = size(Sys.gStrain);

if n1~=nElectrons
  err = sprintf('Sys.gStrain must have %d rows, one per electron spin!',nElectrons);
  return
end

switch n2
  case 1, Sys.gStrain = Sys.gStrain(:,[1 1 1]);
  case 2, Sys.gStrain = Sys.gStrain(:,[1 1 2]);
  case 3 % ok
  otherwise
  err = sprintf('Sys.gStrain must have 1, 2, or 3 columns!');
end

if any(Sys.gStrain(:)<0)
  err = 'Sys.gStrain must contain nonnegative values!';
  return
end


% D and E strain (Sys.DStrain, Sys.DStrainCorr)
%-------------------------------------------------------------------------------
if ~isfield(Sys,'DStrain') || isempty(Sys.DStrain)
  Sys.DStrain = zeros(nElectrons,3);
end
if ~isfield(Sys,'DStrainCorr')
  Sys.DStrainCorr = zeros(1,nElectrons);
end

[n1,n2] = size(Sys.DStrain);

if n1~=nElectrons
  err = sprintf('Sys.DStrain must have %d rows, one per electron spin!',nElectrons);
  return
end

switch n2
  case 1, Sys.DStrain = [Sys.DStrain zeros(nElectrons,2)];
  case 2, Sys.DStrain = [Sys.DStrain zeros(nElectrons,1)];
  case 3 % ok
  otherwise
  err = sprintf('Sys.DStrain must have 1, 2, or 3 columns!');
end

if numel(Sys.DStrainCorr)~=nElectrons
  err = sprintf('Sys.DStrainCorr must contain %d elements, since you have %d electron spins',...
    nElectrons,nElectrons);
end

if any(Sys.DStrain(:)<0)
  err = 'Sys.DStrain must contain nonnegative values!';
  return
end

if any(Sys.DStrainCorr<-1) || any(Sys.DStrainCorr>1)
  err = 'D-E strain correlation coefficient in Sys.DStrainCorr must be between -1 and +1.';
  return
end


% Sys.HStrain, Sys.AStrain
%-------------------------------------------------------------------------------
BroadeningType = {'HStrain','AStrain'};
Elements = [3,3,2];
for k = 1:numel(BroadeningType)
  fld = BroadeningType{k};
  if ~isfield(Sys,fld) || isempty(Sys.(fld))
    Sys.(fld) = zeros(1,Elements(k));
    continue
  end
  p = Sys.(fld);
  if Elements(k)==3
    if numel(p)==1, p = p([1 1 1]); end
    if numel(p)==2, p = p([1 1 2]); end
  end
  if numel(p)~=Elements(k)
    err = sprintf('Sys.%s must have %d elements!',fld,Elements(k));
    if ~isempty(err), return; end
  end
  if any(p(:)<0)
    err = sprintf('Sys.%s must contain nonnegative values!',fld);
    if ~isempty(err), return; end
  end
  Sys.(fld) = p;    
end


% Population vector (Sys.Pop)
%===============================================================================
if ~isfield(Sys,'Pop')
  Sys.Pop = [];
end
if ~isempty(Sys.Pop)
  if ~isvector(Sys.Pop)
    err = 'Sys.Pop must be a row or column vector.';
    return
  end
end


% Diffusion tensor (Sys.tcorr, Sys.logtcorr, Sys.Diff, Sys.logDiff, Sys.DiffFrame)
%===============================================================================
% Rotational correlation time, rotational diffusion rate
fields = {'Diff','logDiff','tcorr','logtcorr'};
for k = 1:numel(fields)
  if isfield(Sys,fields{k})
    if ~any(numel(Sys.(fields{k}))==[0 1 2 3])
      err = sprintf('Sys.%s must have 1, 2, or 3 elements.',fields{k});
      return
    end
  end
end
% Euler angles for diffusion tensor
if isfield(Sys,'DiffFrame')
  err = sizecheck(Sys,'DiffFrame',[1 3]);
  if ~isempty(err); return; end
end


% Multiple Order Zeeman Hamiltonian
%===============================================================================
Sys.MO_present = false;
if any(strncmp('Ham',fieldnames(Sys),3))
  Hamstr = cell(9,9,17);
  for lB = 8:-1:0
    for lS = 8:-1:0
      lmin = abs(lB-lS);
      for l = (lB+lS):-1:lmin
        Hamstr{lB+1,lS+1,(l-lmin)+1} = sprintf('Ham%i%i%i',lB,lS,l);
      end
    end
    Bstr{lB+1} = ['B',num2str(lB)];
  end
  field = isfield(Sys,Hamstr);
  if any(field(:))
    Sys.MO_present = true;
    
    % check for D and Bk
    lB0 = field(1,:,:);
    if any(lB0(:))
      if D_present && field(1,3,1)
        err = 'Cannot use Sys.D and Sys.Ham022 simultaneously. Remove one of them.';
        return
      end
      if any(squeeze(field(1,:,1)).*isfield(Sys,Bstr))
        err = 'Cannot use higher order operators and corresponding general parameters simultaneously. Remove one of them.';
      end
    end
    %check for g
    lB1 = field(2,2,:);
    if any(lB1(:)) && any(Sys.g(:))
      err = 'Cannot use Sys.g and and Sys.Ham112 or Sys.Ham110 simultaneously. Remove one of them.';
    end
    ls =find(field);
    % get l
    [rowsub, colsub, pagsub] = ind2sub([9,9,17], ls);
    l = pagsub - 1 + abs(rowsub-colsub);
    for n = 1:length(ls)
      str = Hamstr{ls(n)};
      if issize(Sys.(str),[nElectrons,1])
        Sys.(str) = [zeros(nElectrons,l(n)), Sys.(str),zeros(nElectrons,l(n))];
      else
        if ~issize(Sys.(str),[nElectrons,2*l(n)+1])
          if ~issize(Sys.(str),[2*l(n)+1,1]) || nElectrons~=1
            err = sprintf('Sys.%s has wrong size!',str);
            return
          else
            Sys.(str) = Sys.(str).';
          end
        end
      end
    end
  end
end


% Orbital angular momentum (Sys.L, Sys.soc, Sys.orf, Sys.CF*)
%===============================================================================
if isfield(Sys,'L') && ~isempty(Sys.L)
  % Guard against invalid type
  if any(~isreal(Sys.L)) || any(mod(real(Sys.L),1)) || any(Sys.S<0)
    err = 'Orbital angular momentum in Sys.L must be nonnegative integers.';
    return
  end
  if numel(Sys.L)~=nElectrons
    err = 'Sys.L and Sys.S must have the same number of elements.';
    return
  end
  if ~isfield(Sys,'soc')
    err = 'Sys.L is given, but no spin-orbit coupling is defined in Sys.soc.';
    return
  end
  if isempty(Sys.soc) || any(~isreal(Sys.soc))
    err = 'Spin-orbit coupling in soc must be real numbers.';
    return
  end
  if size(Sys.soc,1)~=nElectrons
    if issize(Sys.soc,[1,nElectrons])
      Sys.soc = Sys.soc.';
    else
      err = 'Number of spin-orbit couplings in Sys.soc must match number of spins in Sys.S.';
      return
    end
  end
  if ~isfield(Sys,'orf')
    Sys.orf = ones(nElectrons,1);
  else
    if length(Sys.orf)~=nElectrons
      err ='Number of orbital reduction factors must match number of orbital angular momenta!';
      return
    end
    if isempty(Sys.orf) || any(~isreal(Sys.orf))
      err = 'Orbital reduction factors in Sys.orf must be real numbers.';
      return
    end
  end
  for k = 1:12
    fieldname = sprintf('CF%d',k);
    if ~isfield(Sys,fieldname), continue; end
    CFk = Sys.(fieldname);
        
    if size(CFk,1)~=nElectrons
      sn = num2str(nElectrons);
      err = ['Field Sys.', fieldname, ' has to have ',sn,...
        ' rows, since there are ', sn,' orbital angular momenta.'];
    end
    
    if size(CFk,2)==1
      CFk = [zeros(nElectrons,k) CFk(:) zeros(nElectrons,k)];
    elseif size(CFk,2)==2*k+1
      % full form
    else
      err = ['Field Sys.', fieldname, ' has ', num2str(size(CFk,2)), ...
        ' instead of ', num2str(2*k+1),' coloumns.'];
    end
    Sys.(fieldname) = CFk;
  end
else
  if isfield(Sys,'orf') && ~isempty(Sys.orf)
    err = 'Sys.orf is given, but Sys.L is missing. Specify Sys.L.';
    return
  end
  if isfield(Sys,'soc') && ~isempty(Sys.soc)
    err = 'Sys.soc is given, but Sys.L is missing. Specify Sys.L.';
    return
  end
  for k = 1:12
    fn = sprintf('CF%d',k);
    if isfield(Sys,fn) && ~isempty(Sys.(fn))
      err = sprintf('Sys.%s is given, but Sys.L is missing. Specify Sys.L.',fn);
      return
    end
  end
  Sys.L = [];
  Sys.orf = [];
end
Sys.nL = numel(Sys.L);


% Final tasks
%-------------------------------------------------------------------------------
Sys.Spins = [Sys.S(:); Sys.I(:); Sys.L(:)].';
Sys.nStates = hsdim(Sys.Spins);

FullSys = Sys;
FullSys.processed = true;

return
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


%-------------------------------------------------------------------------------
function ok = issize(A,siz)
ok = all(size(A)==siz);
return

%-------------------------------------------------------------------------------
function msg = sizecheck(Sys,FieldName,siz)
ok = all(size(Sys.(FieldName))==siz);
if ok
  msg = '';
else
  msg = sprintf('Spin system field %s has wrong size for the given spins.',FieldName);
end
return

%-------------------------------------------------------------------------------
function err = pa_obsolete_message(Sys,pa,Frame)
if isfield(Sys,pa)
  err = sprintf('Obsolete field Sys.%s. Use Sys.%s instead.',pa,Frame);
else
  err = '';
end
return
