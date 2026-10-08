% strains_parse  Parse and check strain parameter references
%
%   refs = strains_parse(StrainPars,Sys,nElectrons,nNuclei)
%
%   Parses the list of strain parameter references in StrainPars (e.g.
%   {'g(1)','A(1,3)','DFrame(2)'}) against the fields of the raw spin system
%   structure Sys (as entered by the user, before validatespinsys).
%
%   Output: refs, a struct array with one element per reference, with fields
%     Name         reference as given by the user
%     FieldName    name of referenced spin system field
%     Index        linear index into the field
%     Subscripts   [row col] subscripts into the field
%     FieldSize    size of the field
%     Form         input form of the field ('iso','iso1','axial','pv','sym',
%                  'full','D','DE','eeqQ','eeqQeta','dip1','dip2','dip3',
%                  'dvec','J','frame')
%
%   Errors are thrown directly.

function refs = strains_parse(StrainPars,Sys,nElectrons,nNuclei)

% Normalize input to a cell array of char vectors
if isstring(StrainPars) || ischar(StrainPars)
  StrainPars = cellstr(StrainPars);
elseif iscell(StrainPars)
  StrainPars = cellfun(@convertStringsToChars,StrainPars,'UniformOutput',false);
end
if ~iscell(StrainPars) || ~all(cellfun(@(x)ischar(x)&&isrow(x),StrainPars))
  error('Sys.StrainPars must be a cell array of strings, e.g. {''g(1)'',''g(3)''}.');
end

tensorFields = {'g','A','D','Q','ee','J','dip','dvec','nn','sigma'};
frameFields = {'gFrame','AFrame','DFrame','QFrame','eeFrame','nnFrame','sigmaFrame'};

nPars = numel(StrainPars);
refs = struct('Name',{},'FieldName',{},'Index',{},'Subscripts',{},'FieldSize',{},'Form',{});
for p = 1:nPars
  str = strtrim(StrainPars{p});
  % Split into field name and index part, e.g. 'A(1,3)' -> 'A', [1 3]
  fieldName = regexp(str,'^[A-Za-z]\w*','match','once');
  indexStr = strtrim(str(numel(fieldName)+1:end));
  subs = [];
  if ~isempty(indexStr)
    tok = regexp(indexStr,'^\(\s*(\d+)\s*,\s*(\d+)\s*\)$','tokens','once'); % (row,col)
    if isempty(tok)
      tok = regexp(indexStr,'^\(\s*(\d+)\s*\)$','tokens','once'); % (linear index)
    end
    subs = str2double(tok);
  end
  if isempty(fieldName) || (~isempty(indexStr) && isempty(subs))
    error('Sys.StrainPars: ''%s'' is not a valid reference. Use e.g. ''g(1)'' or ''A(1,3)''.',str);
  end
  if ~any(strcmp(fieldName,[tensorFields frameFields]))
    error('Sys.StrainPars: strains for Sys.%s are not supported.',fieldName);
  end
  if ~isfield(Sys,fieldName) || isempty(Sys.(fieldName))
    error('Sys.StrainPars: ''%s'' references Sys.%s, which is not given.',str,fieldName);
  end
  value = Sys.(fieldName);
  siz = size(value);

  % Resolve index
  if isempty(subs)
    if numel(value)~=1
      error('Sys.StrainPars: ''%s'' refers to Sys.%s, which has more than one element. Use an index.',str,fieldName);
    end
    idx = 1;
  elseif isscalar(subs)
    idx = subs;
    if idx<1 || idx>numel(value)
      error('Sys.StrainPars: index in ''%s'' is out of range.',str);
    end
  else
    r = subs(1);
    c = subs(2);
    if r<1 || r>siz(1) || c<1 || c>siz(2)
      error('Sys.StrainPars: index in ''%s'' is out of range.',str);
    end
    idx = sub2ind(siz,r,c);
  end
  [r,c] = ind2sub(siz,idx);

  % Determine input form of the field
  isFrame = any(strcmp(fieldName,frameFields));
  if isFrame
    tensorName = fieldName(1:end-5);
    tensorForm = inputform(Sys,tensorName,nElectrons,nNuclei);
    if any(strcmp(tensorForm,{'full','sym'}))
      error('Sys.StrainPars: ''%s'' is not allowed, since Sys.%s is given as a full matrix.',str,tensorName);
    end
    form = 'frame';
  else
    form = inputform(Sys,fieldName,nElectrons,nNuclei);
    if strcmp(form,'full')
      error('Sys.StrainPars: ''%s'' is not allowed, since Sys.%s is given as a full 3x3 matrix. Use the [xx yy zz xy xz yz] form instead.',str,fieldName);
    end
  end

  % Check for duplicates
  for q = 1:numel(refs)
    if strcmp(refs(q).FieldName,fieldName) && refs(q).Index==idx
      error('Sys.StrainPars: ''%s'' and ''%s'' refer to the same parameter.',refs(q).Name,str);
    end
  end

  refs(p).Name = str;
  refs(p).FieldName = fieldName;
  refs(p).Index = idx;
  refs(p).Subscripts = [r c];
  refs(p).FieldSize = siz;
  refs(p).Form = form;
end

end


%-------------------------------------------------------------------------------
% Determine input form of a tensor field, following the order of size checks
% in validatespinsys.
function form = inputform(Sys,fieldName,nEl,nNuc)

if ~isfield(Sys,fieldName) || isempty(Sys.(fieldName))
  form = '';
  return
end
v = Sys.(fieldName);
siz = size(v);
nElPairs = nEl*(nEl-1)/2;
nNucPairs = nNuc*(nNuc-1)/2;
is = @(s) isequal(siz,s);

switch fieldName
  case 'g'
    form = pvform(v,nEl);
  case 'sigma'
    form = pvform(v,nNuc);
  case 'D'
    if is([3*nEl 3]), form = 'full';
    elseif numel(v)==nEl, form = 'D';
    elseif is([nEl 2]), form = 'DE';
    elseif is([nEl 3]), form = 'pv';
    elseif is([nEl 6]), form = 'sym';
    else, form = 'invalid';
    end
  case 'A'
    if is([3*nNuc 3*nEl]), form = 'full';
    elseif is([1 nNuc]) && nEl==1, form = 'iso1';
    elseif is([nNuc nEl]), form = 'iso';
    elseif is([nNuc 2*nEl]), form = 'axial';
    elseif is([nNuc 3*nEl]), form = 'pv';
    elseif is([nNuc 6*nEl]), form = 'sym';
    else, form = 'invalid';
    end
  case 'Q'
    if is([nNuc 6]), form = 'sym';
    elseif is([3*nNuc 3]), form = 'full';
    elseif numel(v)==nNuc, form = 'eeqQ';
    elseif is([nNuc 2]), form = 'eeqQeta';
    elseif is([nNuc 3]), form = 'pv';
    else, form = 'invalid';
    end
  case {'ee','nn'}
    if strcmp(fieldName,'ee'), nPairs = nElPairs; else, nPairs = nNucPairs; end
    if numel(v)==nPairs, form = 'iso';
    elseif is([nPairs 6]), form = 'sym';
    elseif is([3*nPairs 3]), form = 'full';
    elseif is([nPairs 3]), form = 'pv';
    else, form = 'invalid';
    end
  case 'J'
    form = 'J';
  case 'dip'
    if numel(v)==nElPairs, form = 'dip1';
    else, form = sprintf('dip%d',siz(2));
    end
  case 'dvec'
    form = 'dvec';
  otherwise
    form = 'invalid';
end

end

%-------------------------------------------------------------------------------
function form = pvform(v,n)
siz = size(v);
if isequal(siz,[3*n 3]), form = 'full';
elseif numel(v)==n, form = 'iso';
elseif isequal(siz,[n 2]), form = 'axial';
elseif isequal(siz,[n 3]), form = 'pv';
elseif isequal(siz,[n 6]), form = 'sym';
else, form = 'invalid';
end
end
