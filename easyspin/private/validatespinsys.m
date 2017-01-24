% validatespinsys   Validation of spin system structure
%
%   [FullSys,err] = validatespinsys(Sys)
%
%   Returns an error string in err if spin system Sys is not valid.
%   FullSys is the processed spin system. All missing optional fields
%   are supplemented, and other fields are added:
%
%     nElectrons, nNuclei, Spins, nStates
%     fullg, fullA, fullD
%     I, gn

% Spin system structure fields
%-------------------------------------------------------
% Required fields:
%   S
%   g
%   ee, if more than one electron present
%   A, if nuclei present
% Optional fields (initialized to zero if missing):
%   all others 
%-------------------------------------------------------

function [FullSys,err] = validatespinsys(Sys)

FullSys = [];
err = '';

if ~isstruct(Sys) || (numel(Sys)~=1)
  err = 'Spin system must be a structure.';
  return
end

if isfield(Sys,'processed')
  if Sys.processed
    FullSys = Sys;
    return;
  end
end

%-------- mex compilation check --------------------------------
fileName = 'cubicsolve';
if (exist(fileName,'file')~=3)
 easyspincompile;
  if (exist(fileName,'file')~=3)
    error('EasySpin: Generation of mex files failed.');
  end
end

%-------- spellcheck fields (lower/upper case) ------------------
correctFields = {'S','Nucs','Abund','n',...
  'g','g_','D','ee','ee2','A','A_','Q',...
  'gFrame','DFrame','eeFrame','AFrame','QFrame',...
  'gStrain','HStrain','AStrain','DStrain',...
  'aF','B0','B2','B4','B6','B8','B10','B12',...
  'lw','lwpp','lwEndor','tcorr','logtcorr','Diff','logDiff'};


givenFields = fieldnames(Sys);
for f = 1:numel(givenFields)
  givField = givenFields{f};
  idx = find(strcmpi(givField,correctFields));
  % check if there is a case-insensitive match
  if idx
    corrField = correctFields{idx};
    if ~strcmp(givField,corrField)
      % Wrong capitalization
      error('\n  Fix capitalization: Sys.%s should be Sys.%s\n',givField,corrField);
    else
      % Correct capitalization
    end
  end
end

for ind = find((strncmpi(givenFields,'Ham',3)))
  if isempty(ind), break; end;
  field = givenFields{ind};
  if length(field)~= 6 
    if str2num(field(4))+str2num(field(5))<10  
      error('Wrong length of Sys.%s entry, should be Hamxyz (with x,y,z integer numbers)',field);
    else
      if length(field)~= 7
        error('Wrong length of Sys.%s entry, should be Hamxyz (with x,y,z integer numbers)',field);
      end
    end
  end
  if ~strncmp(field,'Ham',3)
      % Wrong capitalization
      error('Fix capitalization: Sys.%s should be Sys.%s',field,['Ham', field(4:end)]);
  end
end

% -- electron spins field: System.S -------------------------------------
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

% -------- g matrix fields -----------------------------------------
if isfield(Sys,'g')
  
  if isfield(Sys,'g_')
    error('Sys.g and Sys.g_ are given. Remove one of them.');
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
    err = ('Sys.g has wrong size.');
    return
  end

elseif isfield(Sys,'g_')

  Sys.fullg = 0;
  if Sys.fullg
    % full g tensors
  elseif numel(Sys.g_)==nElectrons
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
  %error('Sys.g is missing.');
  if any(strncmp(fieldnames(Sys),'Ham',3))
    Sys.g = 0;
  else
    Sys.g = gfree*ones(nElectrons,3);
  end
  Sys.fullg = 0;
end

% Euler angles for g tensor(s)
if isfield(Sys,'gpa')
  err = sizecheck(Sys,'gpa',[nElectrons 3]);
  disp('***********************************************************************');
  disp('**   Sys.gpa is obsolete. Please use Sys.gFrame instead.             **');
  disp('**   Here is how to convert:                                         **');
  disp('**   If you had                                                      **');
  disp('**      Sys.gpa = [10 -20 56]*pi/180                                 **');
  disp('**   then use                                                        **');
  disp('**      Sys.gFrame = [-56 20 -10]*pi/180                             **');
  disp('**   Change the signs of all three angles, and reverse the order.    **');
  disp('**   For more help, check the documentation on coordinate frames.    **');
  disp('***********************************************************************');
  if ~isempty(err); return; end
  Sys.gFrame = -Sys.gpa(:,[3 2 1]);
end
if ~isfield(Sys,'gFrame') || isempty(Sys.gFrame)
  Sys.gFrame = zeros(nElectrons,3);
end
err = sizecheck(Sys,'gFrame',[nElectrons 3]);
if ~isempty(err); return; end

%---------- zero-field splitting ------------------------------------------

% Supplement partial or missing D
if ~isfield(Sys,'D')
  Sys.D = zeros(nElectrons,3);
  Sys.fullD = 0;
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
if isfield(Sys,'Dpa')
  err = sizecheck(Sys,'Dpa',[nElectrons 3]);
  disp('***********************************************************************');
  disp('**   Sys.Dpa is obsolete. Please use Sys.DFrame instead.             **');
  disp('**   Here is how to convert:                                         **');
  disp('**   If you had                                                      **');
  disp('**      Sys.Dpa = [10 -20 56]*pi/180                                 **');
  disp('**   then use                                                        **');
  disp('**      Sys.DFrame = [-56 20 -10]*pi/180                             **');
  disp('**   Change the signs of all three angles, and reverse the order.    **');
  disp('**   For more help, check the documentation on coordinate frames.    **');
  disp('***********************************************************************');
  if ~isempty(err); return; end
  Sys.DFrame = -Sys.Dpa(:,[3 2 1]);
end
if ~isfield(Sys,'DFrame') || isempty(Sys.DFrame)
  Sys.DFrame = zeros(nElectrons,3);
end
err = sizecheck(Sys,'DFrame',[nElectrons 3]);
if ~isempty(err); return; end


%---------- high-order terms ------------------------------------------

if isfield(Sys,'aF')
  err = sizecheck(Sys,'aF',[1 2]);
  if ~isempty(err), return; end
else
  Sys.aF = [0 0];
end

% B1, B2, B3, B4, B5, B6, etc.
D_present = any(Sys.D(:));
aF_present = any(Sys.aF(:));
for k=1:12
  fieldname = sprintf('B%d',k);
  if ~isfield(Sys,fieldname), continue; end
  Bk = Sys.(fieldname);
  
  if (k==2) && any(Bk(:)) && D_present
    error('Cannot use Sys.D and Sys.B2 simultaneously. Remove one of them.');
  end
  if (k==4) && any(Bk(:)) && aF_present
    error('Cannot use Sys.aF and Sys.B4 simultaneously. Remove one of them.');
  end
  
  if (size(Bk,1)~=nElectrons)
    error('Field Sys.%s has to have %d rows, since there are %d electron spins.',fieldname,nElectrons,nElectrons);
  end
  
  if (size(Bk,2)==1)
    Bk = [zeros(nElectrons,k) Bk(:) zeros(nElectrons,k)];
  elseif (size(Bk,2)==k+1)
    Bk = [Bk zeros(nElectrons,k)];
  elseif (size(Bk,2)==2*k+1)
    % full form
  else
    error('Field Sys.%s has %d instead of %d columns.',fieldname,size(Bk,2),2*k+1);
  end
  
  Sys.(fieldname) = Bk;

end


%---------- electron-electron ------------------------------------------
Sys.fullee = false;
if (nElectrons>1)
  
  eeMatrix = isfield(Sys,'ee');
  JdD = isfield(Sys,'J') || isfield(Sys,'dvec') || isfield(Sys,'eeD');
  
  if ~eeMatrix && ~JdD
    err = 'Spin system contains 2 or more electron spins, but coupling terms are missing (ee; or J, dvec, eeD)!';
    return
  end
  
  if eeMatrix && JdD
    err = 'Both Sys.ee and (Sys.J,Sys.dvec,Sys.eeD) are given - use only one or the other!';
    return
  end
  
  nPairs = nElectrons*(nElectrons-1)/2;
  
  if eeMatrix
    % Bilinear coupling defined via Sys.ee
    %----------------------------------------------------------------------
    
    % Expand isotropic couplings into 3 equal principal values
    if numel(Sys.ee)==nPairs
      Sys.ee = Sys.ee(:)*[1 1 1];
    end
    
    fullee = issize(Sys.ee,[3*nPairs,3]);
    Sys.fullee = fullee;
    if ~fullee
      err = sizecheck(Sys,'ee',[nPairs 3]);
      if ~isempty(err), return; end
    end
    
    if isfield(Sys,'eepa')
      err = sizecheck(Sys,'eepa',[nPairs 3]);
      disp('*********************************************************************');
      disp('**   Sys.eepa is obsolete. Please use Sys.eeFrame instead.         **');
      disp('**   Here is how to convert:                                       **');
      disp('**   If you had                                                    **');
      disp('**      Sys.eepa = [10 -20 56]*pi/180                              **');
      disp('**   then use                                                      **');
      disp('**      Sys.eeFrame = [-56 20 -10]*pi/180                          **');
      disp('**   Change the signs of all three angles, and reverse the order.  **');
      disp('**   For more help, check the documentation on coordinate frames.  **');
      disp('*********************************************************************');
      if ~isempty(err); return; end
      Sys.eeFrame = -Sys.eepa(:,[3 2 1]);
    end
    
  else
    % Bilinear coupling defined via J, dvec and eeD
    %----------------------------------------------------------------------
    % J:    isotropic exchange +J*S1*S2
    % dvec: antisymmetric exchange dvec.(S1xS2)
    % eeD:  dipolar coupling S1.diag(eeD).S2 or S1.eeD.S2
    fullee = false;
    
    % Size check on list of isotropic exchange coupling constants
    if ~isfield(Sys,'J'), Sys.J = zeros(1,nPairs); end
    err = sizecheck(Sys,'J',[1 nPairs]);
    if ~isempty(err), return; end
    
    % Size check on list of antisymmetric exchange vectors
    if ~isfield(Sys,'dvec'), Sys.dvec = zeros(nPairs,3); end
    err = sizecheck(Sys,'dvec',[nPairs,3]);
    if ~isempty(err), return; end
    
    % Size check on dipolar tensor diagonals
    if ~isfield(Sys,'eeD'), Sys.eeD = zeros(nPairs,3); end
    err = sizecheck(Sys,'eeD',[nPairs,3]);
    if ~isempty(err), return; end
    
    % Assert zero traces of dipolar tensors
    if any(sum(Sys.eeD,2))
      err = 'Sys.eeD contains dipolar tensors with non-zero trace. Use Sys.J for this.';
    end
    if ~isempty(err), return; end
    
    % Combine (Sys.J,Sys.dvec,Sys.eeD) into full interaction matrix in Sys.ee
    Sys.fullee = true;
    idx = 1:3;
    for iPair = 1:nPairs
      J = Sys.J(iPair);
      d = Sys.dvec(iPair,:);
      ee = J*eye(3) + ...
         [0 d(3) -d(2); -d(3) 0 d(1); d(2) -d(1) 0] + ...
         diag(Sys.eeD(iPair,:));
      Sys.ee(idx,:) = ee;
      idx = idx + 3;
    end
    
  end
  
  % Check for eeFrame, and supplement or error if necessary
  if fullee
    if isfield(Sys,'eeFrame')
      err = sprintf('Full matrices are specified in ee, so eeFrame is not allowed.');
      if ~isempty(err), return; end
    end
  else
    if ~isfield(Sys,'eeFrame'), Sys.eeFrame = zeros(nPairs,3); end
    err = sizecheck(Sys,'eeFrame',[nPairs 3]);
    if ~isempty(err), return; end
  end

  
  % Isotropic biquadratic exchange term
  %------------------------------------------------------------------------
  if isfield(Sys,'ee2')
    Sys.ee2 = Sys.ee2(:);
    err = sizecheck(Sys,'ee2',[nPairs 1]);
    if ~isempty(err), return; end
  else
    Sys.ee2 = zeros(1,nPairs);
  end
  
end

% Nuclear spins
%==========================================================================

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

if ~isempty(I)
  if any(I==0)
    err = 'System contains nuclei with spin 0. Please remove them.';
    %if ~isempty(err), return; end
  end
end

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
    err = ('Insufficient number of elements in gnscale field of spin system.');
    if ~isempty(err), return; end
  end
else
  Sys.gnscale = ones(1,nNuclei);
end

% ------------------- Hyperfine couplings -------------------------
Sys.fullA = 0;
if (nNuclei>0)
  
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
    elseif issize(Sys.A_,[nNuclei,2*nElectrons]);
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
    
    if issize(Sys.A,[3*nNuclei,3*nElectrons]);
      % Full A matrices
      Sys.fullA = 1;
    elseif issize(Sys.A,[1 nNuclei])
      % Allow simple one-row syntax in the case of 1 eletron spin
      if (nElectrons==1)
        Sys.A = Sys.A.'*[1 1 1];
      else
        err = 'Size of Sys.A matrix is inconsistent with number of electrons and nuclei.';
        if ~isempty(err), return; end
      end
    elseif issize(Sys.A,[nNuclei,nElectrons])
      % Expand isotropic A into 3 equal principal values
      Sys.A = kron(Sys.A,[1 1 1]);
    elseif issize(Sys.A,[nNuclei,2*nElectrons]);
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
  
  if isfield(Sys,'Ascale')
    if numel(Sys.Ascale)<nNuclei
      err = ('Insufficient number of elements in Ascale field of spin system.');
      if ~isempty(err), return; end
    end
  else
    Sys.Ascale = ones(1,nNuclei);
  end

  % Euler angles for A tensor(s)
  if isfield(Sys,'Apa')
    err = sizecheck(Sys,'Apa',[nNuclei,3*nElectrons]);
    disp('*********************************************************************');
    disp('**   Sys.Apa is obsolete. Please use Sys.AFrame instead.           **');
    disp('**   Here is how to convert:                                       **');
    disp('**   If you had                                                    **');
    disp('**      Sys.Apa = [10 -20 56]*pi/180                               **');
    disp('**   then use                                                      **');
    disp('**      Sys.AFrame = [-56 20 -10]*pi/180                           **');
    disp('**   Change the signs of all three angles, and reverse the order.  **');
    disp('**   For more help, check the documentation on coordinate frames.  **');
    disp('*********************************************************************');
    if ~isempty(err); return; end
    idx = reshape(flipud(reshape(idx,3,[])),1,[]);
    Sys.AFrame = -Sys.Apa(:,idx);
  end
  if ~isfield(Sys,'AFrame') || isempty(Sys.AFrame)
    Sys.AFrame = zeros(nNuclei,3*nElectrons);
  end
  err = sizecheck(Sys,'AFrame',[nNuclei,3*nElectrons]);
  if ~isempty(err); return; end
  
end



%------ Nuclear quadrupole ---------------------------
Sys.fullQ = 0;
if (nNuclei>0)
  
  if ~isfield(Sys,'Q')

    Sys.Q = zeros(nNuclei,3);
    Sys.fullQ = 0;

  else

    Sys.fullQ = issize(Sys.Q,[3*nNuclei 3]);
    if ~Sys.fullQ
      if numel(Sys.Q)==nNuclei
        Sys.Q = Sys.Q(:);
        Sys.Q(:,2) = 0;
      end
      if issize(Sys.Q,[nNuclei,2])
        % Convert eeqQ,eta -> principal values
        for iN = 1:nNuclei
          I = Sys.I(iN);
          if (I<1)
            Sys.Q(iN,1:3) = 0;
          else
            eeqQ = Sys.Q(iN,1); % (actually eeqQ/h, in MHz)
            eta = Sys.Q(iN,2);
            Sys.Q(iN,1:3) = eeqQ/(4*I*(2*I-1)) * [-1+eta, -1-eta, 2];
          end
        end
      end
      err = sizecheck(Sys,'Q',[nNuclei,3]);
      if ~isempty(err), return; end
    end
  end

  if isfield(Sys,'Qscale')
    if numel(Sys.Qscale)<nNuclei
      err = ('Insuffient number of elements in Qscale field of spin system.');
      if ~isempty(err), return; end
    end
  else
    Sys.Qscale = ones(1,nNuclei);
  end

  %--------------------
  
  if Sys.fullQ
    for iNuc=1:nNuclei
      Q_ = Sys.Q(3*(iNuc-1)+(1:3),:);
      if norm(Q_-Q_.')/norm(Q_)>1e-5
        err = 'Sys.Q contains asymmetric full Q matrix. Only symmetric Q matrices are allowed.';
        if ~isempty(err), return; end
      end
    end
  end

  % Euler angles for Q tensor(s)
  if isfield(Sys,'Qpa')
    err = sizecheck(Sys,'Qpa',[nNuclei 3]);
    disp('*********************************************************************');
    disp('**   Sys.Qpa is obsolete. Please use Sys.QFrame instead.           **');
    disp('**   Here is how to convert:                                       **');
    disp('**   If you had                                                    **');
    disp('**      Sys.Qpa = [10 -20 56]*pi/180                               **');
    disp('**   then use                                                      **');
    disp('**      Sys.QFrame = [-56 20 -10]*pi/180                           **');
    disp('**   Change the signs of all three angles, and reverse the order.  **');
    disp('**   For more help, check the documentation on coordinate frames.  **');
    disp('*********************************************************************');
    if ~isempty(err); return; end
    Sys.QFrame = -Sys.Qpa(:,[3 2 1]);
  end
  if ~isfield(Sys,'QFrame') || isempty(Sys.QFrame)
    Sys.QFrame = zeros(nNuclei,3);
  end
  err = sizecheck(Sys,'QFrame',[nNuclei 3]);
  if ~isempty(err); return; end

end

% Remove spin-zero nuclei
rmv = Sys.I==0;

if any(rmv)

  if Sys.fullA
    for iNuc = numel(rmv):-1:1
      if ~rmv(iNuc), continue; end
      Sys.A((iNuc-1)*3+(1:3),:) = [];
    end
  else
    Sys.A(rmv,:) = [];
    Sys.AFrame(rmv,:) = [];
  end
    
  Sys.Nucs(rmv) = [];  
  Sys.Q(rmv,:) = [];
  Sys.I(rmv) = [];
  Sys.gn(rmv) = [];
  Sys.QFrame(rmv,:) = [];
  Sys.nNuclei = numel(Sys.gn);
  Sys.Ascale(rmv) = [];
  Sys.Qscale(rmv) = [];
  Sys.gnscale(rmv) = [];
  Sys.n(rmv) = [];
end


% Broadenings (Strains and convolution line widths)
%==========================================================================

% Check Sys.lw
if ~isfield(Sys,'lw'), Sys.lw = [0 0]; end
if numel(Sys.lw)==1, Sys.lw(2) = 0; end
if numel(Sys.lw)~=2, err = ('System.lw has wrong size.'); end
if ~isempty(err), return; end
if any(Sys.lw)<0, err = ('System.lw cannot be negative.'); end
if ~isempty(err), return; end

% Check Sys.lwEndor
if ~isfield(Sys,'lwEndor'), Sys.lwEndor = [0 0]; end
if numel(Sys.lwEndor)==1, Sys.lwEndor(2) = 0; end
if numel(Sys.lwEndor)~=2, err = ('System.lwEndor has wrong size.'); end
if ~isempty(err), return; end
if any(Sys.lwEndor)<0, err = ('System.lwEndor cannot be negative.'); end
if ~isempty(err), return; end

% Check Sys.lwpp
if ~isfield(Sys,'lwpp'), Sys.lwpp = [0 0]; end
if numel(Sys.lwpp)==1, Sys.lwpp(2) = 0; end
if numel(Sys.lwpp)~=2, err = ('System.lwpp has wrong size.'); end
if ~isempty(err), return; end
if any(Sys.lwpp)<0, err = ('Linewidths cannot be negative.'); end
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

% g strain
%------------------------------------------------------------------
if ~isfield(Sys,'gStrain') || isempty(Sys.gStrain)
  Sys.gStrain = zeros(nElectrons,3);
end

[n1,n2] = size(Sys.gStrain);

if (n1~=nElectrons)
  err = sprintf('Sys.gStrain must have %d rows, one per electron spin!',nElectrons);
  return
end

switch n2
  case 1, Sys.gStrain = Sys.gStrain(:,[1 1 1]);
  case 2, Sys.gStrain = Sys.gStrain(:,[1 1 2]);
  case 3, % ok
  otherwise
  err = sprintf('Sys.gStrain must have 1, 2, or 3 columns!');
end

if any(Sys.gStrain(:,1)<0) || any(Sys.gStrain(:,2)<0)
  err = 'Sys.gStrain must contain nonnegative values!';
  return
end

% D and E strain
%------------------------------------------------------------------
if ~isfield(Sys,'DStrain') || isempty(Sys.DStrain)
  Sys.DStrain = zeros(nElectrons,3);
end
if ~isfield(Sys,'DStrainCorr')
  Sys.DStrainCorr = zeros(1,nElectrons);
end

[n1,n2] = size(Sys.DStrain);

if (n1~=nElectrons)
  err = sprintf('Sys.DStrain must have %d rows, one per electron spin!',nElectrons);
  return
end

switch n2
  case 1, Sys.DStrain = [Sys.DStrain zeros(nElectrons,2)];
  case 2, Sys.DStrain = [Sys.DStrain zeros(nElectrons,1)];
  case 3, % ok
  otherwise
  err = sprintf('Sys.DStrain must have 1, 2, or 3 columns!');
end

if numel(Sys.DStrainCorr)~=nElectrons
  err = sprintf('Sys.DStrainCorr must contain %d elements, since you have %d electron spins',...
    nElectrons,nElectrons);
end

if any(Sys.DStrain(:,1)<0) || any(Sys.DStrain(:,2)<0)
  err = 'Sys.DStrain must contain nonnegative values!';
  return
end

if any(Sys.DStrainCorr<-1) || any(Sys.DStrainCorr>1)
  err = 'D-E strain correlation coefficient in Sys.DStrainCorr must be between -1 and +1.';
  return
end


% Check width parameters for correct size, assure that they are not negative
BroadeningType = {'HStrain','AStrain'};
Elements = [3,3,2];
for k = 1:numel(BroadeningType)
  fld = BroadeningType{k};
  if ~isfield(Sys,fld) || isempty(Sys.(fld))
    Sys.(fld) = zeros(1,Elements(k));
    continue
  end
  p = Sys.(fld);
  if (Elements(k)==3)
    if (numel(p)==1), p = p([1 1 1]); end
    if (numel(p)==2), p = p([1 1 2]); end
  end
  if numel(p)~=Elements(k)
    err = sprintf('Sys.%s must have %d elements!',fld,Elements(k));
    if ~isempty(err), return; end
  end
  if any(p<0)
    err = sprintf('Sys.%s must contain nonnegative values!',fld);
    if ~isempty(err), return; end
  end
  Sys.(fld) = p;    
end

%if any(Sys.gStrain) & any(Sys.HStrain)
%  err = sprintf('Sys.gStrain and Sys.HStrain affect the spectrum in the same way. Use only one of the two.');
%end

% Diffusion tensor
%=========================================================================================
% Euler angles for diffusion tensor
if isfield(Sys,'Diffpa')
  err = sizecheck(Sys,'Diffpa',[1 3]);
  disp('***********************************************************************');
  disp('**   Sys.Diffpa is obsolete. Please use Sys.DiffFrame instead.       **');
  disp('**   Here is how to convert:                                         **');
  disp('**   If you had                                                      **');
  disp('**      Sys.Diffpa = [10 -20 56]*pi/180                              **');
  disp('**   then use                                                        **');
  disp('**      Sys.DiffFrame = [-56 20 -10]*pi/180                          **');
  disp('**   Change the signs of all three angles, and reverse the order.    **');
  disp('**   For more help, check the documentation on coordinate frames.    **');
  disp('***********************************************************************');
  if ~isempty(err); return; end
  Sys.DiffFrame = -Sys.Diffpa(:,[3 2 1]);
end

% Multiple Order Hamiltonian
for lB = 0:8
  for lS = 0:8
    for l=abs(lB-lS):(lB+lS)
      str = ['Ham',num2str([lB,lS,l],'%i%i%i')];
      if isfield(Sys,str)
        if lB == 0
          % check for D, aF, and Bk
          if lS == 2 && D_present
            error('Cannot use Sys.D and Sys.Ham022 simultaneously. Remove one of them.');
          end
          if lS == 4 && aF_present
            error('Cannot use Sys.aF and Sys.Ham044 simultaneously. Remove one of them.');
          end          
          Bstr = ['B',num2str(lS)];
          if isfield(Sys,Bstr)
             error('Cannot use Sys.%s and Sys.%s simultaneously. Remove one of them.',Bstr,str);
          end
        end
        if lB == 1 && any(Sys.g(:))
             error('Cannot use Sys.g and Sys.%s simultaneously. Remove one of them.',str);
        end  
        if issize(Sys.(str),[nElectrons,1])
          Sys.(str) = [zeros(nElectrons,l), Sys.(str),zeros(nElectrons,l)];
        else
          if ~issize(Sys.(str),[nElectrons,2*l+1])
            if ~issize(Sys.(str),[2*l+1,1])
              error('Sys.%s has wrong size!',str);
            else
              Sys.(str) = Sys.(str).';
            end
          end
        end
      end
    end
  end
end

  



%--------------------------------------------------------------------------
Sys.Spins = [Sys.S(:); Sys.I(:)].';
Sys.nStates = hsdim(Sys.Spins);

FullSys = Sys;
FullSys.processed = 1;

return
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%------------------
function ok = issize(A,siz)
ok = all(size(A)==siz);
return

%------------------
function msg = sizecheck(Sys,FieldName,siz)
ok = all(size(Sys.(FieldName))==siz);
if ok
  msg = '';
else
  %msg = sprintf('Spin system field %s must be a %dx%d array.',...
  %  FieldName,siz(1),siz(2));
  msg = sprintf('Spin system field %s has wrong size for the given spins.',FieldName);
end
return
