function funcs = esfit_argsparams()

funcs.getparaminfo = @getparaminfo;
funcs.setparamvalues = @setparamvalues;
funcs.getparamvalues = @getparamvalues;
funcs.checkparcompatibility = @checkparcompatibility;
funcs.parinfoequal = @parinfoequal;
funcs.validargs = @validargs;

end

%===============================================================================
% function for parsing vary, lower, or upper-bound inputs
function parinfo = getparaminfo(args)

if ~iscell(args)
  error('Model function arguments are expected as a cell array.')
end
nArgs = numel(args);

p = 1; % parameter index
parinfo = struct;
for a = 1:nArgs
  par = args{a};
  if isempty(par), continue; end
  noCell = ~iscell(par);
  if noCell, par = {par}; end
  for c = 1:numel(par)
    if ~isstruct(par{c})
      error('Arg %d, cell %d: struct expected.',a,c);
    end
    allFields = fieldnames(par{c});
    allValues = struct2cell(par{c});
    for iField = 1:numel(allFields)
      FieldValue = allValues{iField};
      noParCell = ~iscell(FieldValue);
      if noParCell, FieldValue = {FieldValue}; end
      for pc = 1:numel(FieldValue)
        value = FieldValue{pc};
        siz = size(value);
        if ~isfloat(value)
          continue
        end
        idx = 1:numel(value);
        [irow,icol] = ind2sub(size(value),idx);
        idxordered = reshape(idx,size(value)).'; % more intuitive order in output/GUI for array inputs
        for i = idxordered(:).'
          parinfo(p).Arg = a;
          if noCell
            parinfo(p).Cell = 0;
          else
            parinfo(p).Cell = c;
          end
          if noParCell
            parinfo(p).ParCellIndex = 0;
          else
            parinfo(p).ParCellIndex = pc;
          end
          parinfo(p).FieldName = allFields{iField};
          parinfo(p).Index = idx(i); % linear index
          parinfo(p).Subscripts = [irow(i) icol(i)]; % 2D subscripts
          parinfo(p).FieldSize = siz;

          % Build parameter name
          idxName = '';
          if ~isscalar(value)
            idxName = sprintf('(%d,%d)',irow(i),icol(i));
          end
          cellName = '';
          if ~noCell
            cellName = sprintf('{%d}',c);
          end
          ParCellName = '';
          if ~noParCell
            ParCellName = sprintf('{%d}',pc);
          end
          Name = sprintf('arg%1d%s.%s%s%s',a,cellName,allFields{iField},ParCellName,idxName);
          parinfo(p).Name = Name;

          p = p + 1;

        end % index
      end % cell index
    end % field
  end % cell
end % argument

end


%===============================================================================
function args = setparamvalues(args,parinfo,newvals)

nParams = numel(parinfo);
if numel(newvals)~=nParams
    error('Number of elements in parameter vector doesn''t match number of free parameters.');
end

for p = 1:nParams
    a = parinfo(p).Arg;
    c = parinfo(p).Cell;
    pc = parinfo(p).ParCellIndex;
    f = parinfo(p).FieldName;
    i = parinfo(p).Index;
    if c>0
      if pc>0
        args{a}{c}.(f){pc}(i) = newvals(p);
      else
        args{a}{c}.(f)(i) = newvals(p);
      end
    else
      if pc>0
        args{a}.(f){pc}(i) = newvals(p);
      else
        args{a}.(f)(i) = newvals(p);
      end
    end
end

end


%===============================================================================
function vals = getparamvalues(args,parinfo)

nParams = numel(parinfo);
vals = NaN(nParams,1);

for p = 1:nParams
    a = parinfo(p).Arg;
    c = parinfo(p).Cell;
    pc = parinfo(p).ParCellIndex;
    f = parinfo(p).FieldName;
    i = parinfo(p).Index;
    if c>0
      if pc>0
        vals(p) = args{a}{c}.(f){pc}(i);
      else
        vals(p) = args{a}{c}.(f)(i);
      end
    else
      if pc>0
        vals(p) = args{a}.(f){pc}(i);
      else
        vals(p) = args{a}.(f)(i);
      end
    end
end

end


%===============================================================================
function checkparcompatibility(parinfo,args)

nParams = numel(parinfo);
nArgs = numel(args);

for p = 1:nParams
    a = parinfo(p).Arg;
    c = parinfo(p).Cell;
    pc = parinfo(p).ParCellIndex;
    f = parinfo(p).FieldName;
    i = parinfo(p).Index;
    
    if a>nArgs
        error('Parameter %d: Only %d args given, cannot access arg no %d.',p,nArgs,a);
    end
    arg_ = args{a};
    if iscell(arg_)
        if c>numel(arg_)
            error('Parameter %d: Not enough cell entries in argument %d.',p,a);
        end
        arg_ = arg_{c};
    end
    if ~isstruct(arg_)
        error('Structure expected.');
    end
    if ~isfield(arg_,f)
        error('Parameter %d: Field .%s is missing in input %d.',p,f,a);
    end
    if iscell(arg_.(f))
      N = numel(arg_.(f){pc});
    else
      N = numel(arg_.(f));
    end
    if i>N
        error('Parameter %d: Field .%s doesn''t have enough elements - at least %d expected.',p,f,i);
    end
    
end

end

%===============================================================================
function err = parinfoequal(p1,p2)
nParams = numel(p1);
if numel(p2)~=nParams
    error('Number of parameters is not equal.');
end
err = '';
for p = 1:nParams
    if p1(p).Arg~=p2(p).Arg
        err = sprintf('Parameter %d: arg is different (%d vs %d).',p,p1(p).Arg,p2(p).Arg);
    elseif p1(p).Cell~=p2(p).Cell
        err = sprintf('Parameter %d: cell is different (%d vs %d).',p,p1(p).Cell,p2(p).Cell);
    elseif ~strcmp(p1(p).FieldName,p2(p).FieldName)
        err = sprintf('Parameter %d: field name is different (''%s'' vs ''%s'').',p,p1(p).FieldName,p2(p).FieldName);
    elseif p1(p).Index~=p2(p).Index
        err = sprintf('Parameter %d: index is different (%d vs %d).',p,p1(p).Index,p2(p).Index);
    end
    if err, return; end
end
end

%===============================================================================
function validargs(args)

if iscell(args)
  for a = 1:numel(args)
    args_ = args{a};
    if ~iscell(args_) && ~isstruct(args_)
      error('args{%d} must be a cell array or structure array.',a);
    end
    if ~iscell(args_)
      args_ = {args_};
    end
    for c = 1:numel(args_)
      if ~isstruct(args_{c})
        error('args{%d}{%d} must be a structure.',a,c);
      end
    end
  end
elseif isvector(args)
  if min(size(args))>1
    error('Parameter vector must be a row or column array.');
  end
else
  error('args must be a one-dimensional cell array or a numeric array.');
end

end
