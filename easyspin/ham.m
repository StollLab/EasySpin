% ham  Spin Hamiltonian 
%
%   [H0,mux,muy,muz] = ham(Sys)
%   [H0,muzL] = ham(Sys,B0)
%   H = ham(Sys,B0)
%   H = ham(Sys)
%   __ = ham(__,'sparse')
%
%   Constructs the spin Hamiltonian matrix, or its field-independent and
%   field-dependent components, for the spin system Sys.
%
%   Input:
%     Sys       spin system specification structure
%     B0        vector specifying the static magnetic field (in mT) in the
%               molecular frame (e.g. [350 0 0] is along the molecular x axis)
%     'sparse'  if given, sparse instead of full matrices are returned
%
%   Output:
%     H            complete spin Hamiltonian (MHz)
%     H0           field-independent part of spin Hamiltonian (MHz)
%     mux,muy,muz  magnetic dipole moment operators (MHz/mT)
%                  along x, y, z axes of molecular frame
%     muzL         magnetic dipole moment operator (MHz/mT)
%                  along z axis of lab frame
%
%   The spin Hamiltonian components are defined such that
%           H = H0 - B0(1)*mux - B0(2)*muy - B0(3)*muz
%           H = H0 - norm(B0)*muzL

function varargout = ham(SpinSystem,B0,opt)

if nargin==0, help(mfilename); return; end

if nargin<2, B0 = []; end
if nargin<3, opt = ''; end
if ~ischar(opt)
  error('Third argument must be a string, ''sparse'' or ''''.');
end
sparseMatrices = strcmp(opt,'sparse');
if sparseMatrices
  sp = 'sparse';
else
  sp = 'full';
end

[Sys,err] = validatespinsys(SpinSystem);
error(err);

if ~isempty(B0)
  if ~isnumeric(B0) || numel(B0)~=3 || ~isreal(B0)
    error('Magnetic field vector must be 3-element array!');
  end
end

% Build field-independent part of Hamiltonian
H0 = ham_zf(Sys,[],sp) + ...
     ham_ee(Sys,[],sp) + ...
     ham_hf(Sys,[],[],sp) + ...
     ham_nq(Sys,[],sp) + ...
     ham_nn(Sys,[],sp) + ...
     ham_so(Sys,[],sp) + ...
     ham_cf(Sys,[],sp);

% Build magnetic-moment operators
[muxM,muyM,muzM] = ham_ez(Sys,[],sp);  % total electron magnetic moment
if Sys.nNuclei>0
  % Add total nuclear magnetic moment
  [muxMnuc,muyMnuc,muzMnuc] = ham_nz(Sys,[],sp);
  muxM = muxM + muxMnuc;
  muyM = muyM + muyMnuc;
  muzM = muzM + muzMnuc;
end
if Sys.nL>0
  % Add total orbital magnetic moment
  [muxMorb,muyMorb,muzMorb] = ham_oz(Sys,[],sp);
  muxM = muxM + muxMorb;
  muyM = muyM + muyMorb;
  muzM = muzM + muzMorb;
end

% Look for higher-order Zeeman terms
sysfields = fieldnames(Sys);
highest = 0;
higherOrderZeeman = false;
higherOrderZeemanFields = strncmp(sysfields,'Ham',3).';
if any(higherOrderZeemanFields) 
  for n = find(higherOrderZeemanFields)
    if any(Sys.(sysfields{n}))
      higherOrderZeeman = true;
      order = str2double(sysfields{n}(4));
      if order>highest && order<4, highest = order; end
    end
  end
end

if isempty(B0)
  
  if higherOrderZeeman
    % full tensors up to the highest used order will be provided
    zHo = ham_ezho(Sys,[],opt);
    switch highest
      case 0
        if nargout==4
          varargout = {H0+zHo,muxM,muyM,muzM};
        elseif nargout==1
          varargout{1} = {H0+zHo,muxM,muyM,muzM};
        else
          error('1 or 4 output arguments expected!');
        end
      case 1
        if nargout==4
          varargout = {H0+zHo{1},muxM+zHo{2}{1},muyM+zHo{2}{2},muzM+zHo{2}{3}};
        elseif nargout==1
          varargout{1} = {H0+zHo{1},muxM+zHo{2}{1},muyM+zHo{2}{2},muzM+zHo{2}{3}};
        else
          error('1 or 4 output arguments expected!');
        end
      otherwise
        zHo{1} = zHo{1}+H0;
        Gn = {muxM,muyM,muzM};
        for k = 3:-1:1
          zHo{2}{k} = zHo{2}{k}+Gn{k};
        end
        if nargout==highest+1
          varargout = zHo;
        elseif nargout==1
          varargout{1} = zHo;
        else
          error('Wrong number of output arguments!');
        end
    end

  else
    
    if nargout==4
      varargout = {H0,muxM,muyM,muzM};
    elseif nargout==1
      varargout{1} = {{H0,muxM,muyM,muzM}};
    else
      error('1 or 4 output arguments expected!');
    end
  end
  
  return
  
end


% Magnetic field is provided
%------------------------------------------------------------------------
if norm(B0)>0
  nB0 = B0/norm(B0);
else
  nB0 = [0 0 0];
end
muzL = nB0(1)*muxM + nB0(2)*muyM + nB0(3)*muzM;

if higherOrderZeeman
  if nargout==1 || nargout==0
    varargout = {H0 - norm(B0)*muzL + ham_ezho(Sys,B0,[],opt)};
  elseif nargout==2 && highest<2
    zHo = ham_ezho(Sys,[],opt);
    if highest==1
      muzL = muzL + nB0(1)*zHo{2}{1} + nB0(2)*zHo{2}{2} + nB0(3)*zHo{2}{3};
      H0 = H0 + zHo{1};
    else
      H0 = H0 + zHo;
    end
    varargout = {H0,muzL};
  else
    error('Wrong number of output arguments!');
  end
else
  if nargout==1 || nargout==0
    varargout = {H0 - norm(B0)*muzL};
  elseif nargout==2
    varargout = {H0,muzL};
  else
    error('Wrong number of output arguments!');
  end
end

end
