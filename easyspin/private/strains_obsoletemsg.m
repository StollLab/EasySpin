% strains_obsoletemsg  Error message for obsolete strain fields
%
%   msg = strains_obsoletemsg(Sys)
%
%   Returns an error message with replacement code if the raw spin system Sys
%   contains any of the obsolete fields gStrain, AStrain, gAStrainCorr,
%   DStrain, DStrainCorr. Returns '' otherwise.

function msg = strains_obsoletemsg(Sys)

msg = '';
oldFields = {'gStrain','AStrain','gAStrainCorr','DStrain','DStrainCorr'};
given = isfield(Sys,oldFields);
if ~any(given), return; end

if isfield(Sys,'S'), nEl = numel(Sys.S); else, nEl = 1; end
getf = @(f) getfieldvalue(Sys,f);

pars = {};   % strain parameter references
fwhm = [];   % widths
code = {};   % replacement code lines for spin system fields
gIdx = zeros(1,3); % position of g(1,j) in pars, for g-A correlation
AIdx = zeros(1,3); % position of A(1,j) in pars

% g strain: independent widths along the principal axes of each g tensor
gs = getf('gStrain');
if any(gs(:))
  gs = expand3(gs,nEl);
  g = getf('g');
  if isempty(g), g = gfree*ones(nEl,1); end
  if ~isequal(size(g),[nEl 3]) || ~isfield(Sys,'g')
    if numel(g)==nEl, g = g(:)*[1 1 1];
    elseif isequal(size(g),[nEl 2]), g = g(:,[1 1 2]);
    else, g = []; % full g: no replacement possible
    end
    if ~isempty(g), code{end+1} = sprintf('Sys.g = %s;',mat2str(g,8)); end
  end
  for e = 1:nEl
    for j = 1:3
      if gs(e,j)~=0
        pars{end+1} = sprintf('g(%d,%d)',e,j); %#ok<*AGROW>
        fwhm(end+1) = gs(e,j);
        if e==1, gIdx(j) = numel(pars); end
      end
    end
  end
end

% A strain: first nucleus, first electron, along the principal axes of A
As = getf('AStrain');
if any(As(:))
  As = expand3(As(:).',1);
  A = getf('A');
  nNuc = 0;
  if isfield(Sys,'Nucs') && ~isempty(Sys.Nucs), nNuc = numel(nucstring2list(Sys.Nucs)); end
  if nEl==1 && nNuc>0 && ~isequal(size(A),[nNuc 3])
    if numel(A)==nNuc, A = A(:)*[1 1 1];
    elseif isequal(size(A),[nNuc 2]), A = A(:,[1 1 2]);
    else, A = [];
    end
    if ~isempty(A), code{end+1} = sprintf('Sys.A = %s;',mat2str(A,8)); end
  end
  for j = 1:3
    if As(j)~=0
      pars{end+1} = sprintf('A(1,%d)',j);
      fwhm(end+1) = As(j);
      AIdx(j) = numel(pars);
    end
  end
end

% Correlation between g and A strain (old default: +1)
corrList = zeros(0,3); % [i j c] for each nonzero correlation coefficient
if any(gIdx) && any(AIdx)
  corr = getf('gAStrainCorr');
  if isempty(corr), corr = 1; end
  for j = 1:3
    if gIdx(j) && AIdx(j)
      corrList(end+1,:) = [gIdx(j) AIdx(j) sign(corr)];
    end
  end
end

% D strain: [FWHM_D FWHM_E] per electron, with D-E correlation
Ds = getf('DStrain');
if any(Ds(:))
  if size(Ds,2)==1, Ds(:,2) = 0; end
  Dcorr = getf('DStrainCorr');
  if isempty(Dcorr), Dcorr = zeros(1,nEl); end
  D = getf('D');
  % D(e,1) and D(e,2) must refer to D and E, so convert other input forms
  isDonly = numel(D)==nEl;
  needDE = ~isequal(size(D),[nEl 2]) && ~(isDonly && ~any(Ds(:,2)));
  if needDE
    if isDonly
      D = [D(:) zeros(nEl,1)];
    elseif isequal(size(D),[nEl 3])
      D = [D(:,3)-(D(:,1)+D(:,2))/2, (D(:,1)-D(:,2))/2];
    else
      D = [];
    end
    if ~isempty(D), code{end+1} = sprintf('Sys.D = %s;',mat2str(D,8)); end
  end
  for e = 1:nEl
    iD = 0;
    if Ds(e,1)~=0
      pars{end+1} = sprintf('D(%d,1)',e); fwhm(end+1) = Ds(e,1); iD = numel(pars);
    end
    if Ds(e,2)~=0
      pars{end+1} = sprintf('D(%d,2)',e); fwhm(end+1) = Ds(e,2);
      if iD && Dcorr(e)~=0
        corrList(end+1,:) = [iD numel(pars) Dcorr(e)];
      end
    end
  end
end

% Correlation matrix
C = eye(numel(pars));
for k = 1:size(corrList,1)
  C(corrList(k,1),corrList(k,2)) = corrList(k,3);
  C(corrList(k,2),corrList(k,1)) = corrList(k,3);
end

% Assemble message
str = sprintf('''%s'',',pars{:});
code{end+1} = sprintf('Sys.StrainPars = {%s};',str(1:end-1));
code{end+1} = sprintf('Sys.StrainFWHM = %s;',mat2str(fwhm,8));
if ~isempty(corrList)
  code{end+1} = sprintf('Sys.StrainCorr = %s;',mat2str(C,8));
end

oldGiven = sprintf('Sys.%s, ',oldFields{given});
msg = sprintf(['%s are no longer supported. Use Sys.StrainPars together with ' ...
  'Sys.StrainFWHM/Sys.StrainCorr or Sys.StrainModes instead.'],oldGiven(1:end-2));
if ~isempty(pars)
  msg = sprintf('%s Replacement:\n  %s',msg,strjoin(code,'\n  '));
end

end

%-------------------------------------------------------------------------------
function v = getfieldvalue(Sys,f)
if isfield(Sys,f), v = Sys.(f); else, v = []; end
end

%-------------------------------------------------------------------------------
% Expand 1, 2, or 3 values per row to 3 principal values
function v = expand3(v,nRows)
if numel(v)==nRows, v = v(:); end
switch size(v,2)
  case 1, v = v(:,[1 1 1]);
  case 2, v = v(:,[1 1 2]);
end
end
