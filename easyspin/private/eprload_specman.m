function [Data,Abscissa,Parameters] = eprload_specman(FileName)

%  eprload_specman Read data from .d01/.exp SpecMan data files
%
%   y = eprload_specman(FileName)
%   [x,y] = eprload_specman(FileName)
%   [x,y,Pars] = eprload_specman(FileName)
%
%   Reads data and axis information from Specman
%   .d01 files.  
%
% Based on original code by Boris Epel & Alexey Silakov


fid=fopen(char(FileName),'r', 'ieee-le');
if fid<1, error(['File ''',FileName,''' can not be open for read.']);end

ndim1=fread(fid, 1,'uint32');       % number of headers, re/im etc.
dformat=fread(fid,1,'uint32');     % format:0-double,1-float
if dformat==1
    sformat='float32';
else
    sformat='double';
end

dstrms = {};
ntotal = 1;

for k=1:ndim1
    ndim2       = fread(fid,1,'int32');
    dstrms{end+1}.dim = fread(fid,4,'int32');
    dstrms{end}.dim(ndim2+1:end) = 1;
    dstrms{end}.first = ntotal;
    dstrms{end}.total = fread(fid,1,'int32');
    ntotal = ntotal + dstrms{end}.total;
end
tmpdat=fread(fid,ntotal,sformat);
fclose(fid);

switch(ndim1)
    case 0,
        error('No data present');
    case 2,
        Data=tmpdat(dstrms{1}.first:(dstrms{1}.first+dstrms{1}.total-1))+...
            1i*tmpdat((dstrms{2}.first:dstrms{2}.first+dstrms{2}.total-1));
        Data=reshape(Data,dstrms{1}.dim');
    case 1,
        Data=reshape(tmpdat,dstrms{1}.dim');
    otherwise,
        %% find if all data have the same dimensions
        dim = dstrms{1}.dim;
        isthesame = 1;
        is1D      = (sum(dim~=1) == 1);
        for k=2:ndim1
            if sum(dim == dstrms{k}.dim) ~= 4,
                isthesame = false;
                break;
            elseif is1D && sum(dstrms{k}.dim~=1)~=1
                is1D = false;
            end
        end
        if isthesame && is1D,
            % read all as columns
            xdim = dim(dim~=1);
            Data=reshape(tmpdat, xdim, ndim1);
        else
            Data=tmpdat;
        end
end

dscname=strrep(FileName, 'd01', 'exp');


Absicca = [];
Parameters = SpecMandsc(dscname);
ax = SpecManpar(Parameters);

if ~isfield(ax, 'x') || size(ax.x, 1)~=size(Data, 1)
    Abscissa{1} = [1:size(Data, 1)];
else
    Abscissa{1} = ax.x;
    Parameters.xlabel = ax.xlabel;
end

if isfield(ax, 'y') 
    if size(ax.y, 1)~=size(Data, 2)
        Abscissa{2} = [1:size(Data, 2)];
    else
        Abscissa{2} = ax.y;
        Parameters.ylabel = ax.ylabel;
    end
end

if (numel(Abscissa)==1)
  Abscissa = Abscissa{1}(:);
end

return

function par = SpecMandsc(filename)
h = fopen(filename);
if h<0
    disp('Description file was not found.');
    par = [];
    return;
end
olda = '';
par = [];
forbidden = '~!@#$%^&*()./\';
section = '';
text = 1; prg = 1;
while feof(h)<1
  s = strtrim(fgetl(h));
  if isempty(s), continue; end;
  sect = find(s=='[' | s==']');
  
  %   this is a section header
  if size(sect, 2)==2 && sect(1)==1
    section = s(sect(1)+1:sect(2)-1);
    section(section=='-')='';
  else
    switch section
      case 'text'
        par.(['text', num2str(text)]) = s;
        text = text + 1;
      case 'program'
        par.(['prg', num2str(prg)]) = s;
        prg = prg + 1;
      otherwise
        [a,s]=strtok(s, '=');
        a = strtrim(a);
        a(a=='/' | a=='\' | a==' ')='_';
        par.([section,'_',a]) = s(2:end);
    end
  end
end
fclose(h);
return

function res = SpecManpar(par)

prefix = ['n', 'u', 'm', 'k', 'M', 'G'];
koeff  = [1E-9, 1E-6, 1E-3, 1E3, 1E6, 1E9];
res.title = safeget(par, 'general_name', '?');
sweepax = {};
key = 'transient';
fullfield = ['sweep_', key];
idx = 0;
triggers = str2num(safeget(par, 'streams_triggers', '1'));
while isfield(par, fullfield)
    [ax.t, str] = strtok(getfield(par, fullfield), ',');
    ax.t = strtrim(ax.t);
    [ax.size, str] = strtok(str(2:end), ',');
    if ax.t=='S' | ax.t=='I' | ax.t=='A' | ax.t=='R', ax.size = 1; else ax.size = str2num(ax.size); end
    [ax.reps, str] = strtok(str(2:end), ',');
    ax.reps = str2num(ax.reps);
    ax.var = {};
    while ~isempty(str)
        [ax.var{end+1}, str] = strtok(str(2:end), ',');
    end
    sweepax{end+1,1} = ax;
    fullfield = ['sweep_sweep', num2str(idx)];
    idx = idx +1;
end
sweepax{1}.size=sweepax{1}.size*triggers;
res.sweepax = sweepax;

axislabel = 'xyz';
counter = 1;
for k = 1:size(sweepax, 1)
    arr = [];
    asize = sweepax{k}.size;
    if asize > 1
        switch sweepax{k}.t
            case {'I', 'A'}
                tempparam = 'trans';
                par.params_trans = '1sl step 1sl;';
            case 'T'
                tempparam = 'trans';
                dwell_time_str = safeget(par, 'streams_dwelltime', '1 ns');
                dwell_time_str = strtrim(gettoken(dwell_time_str, ','));
                par.params_trans = ['0 ns step ', dwell_time_str,';'];
            otherwise
                tempparam = sweepax{k}.var{1};
                tempparam(strfind(tempparam, ' ')) = '_';
        end
        % check if this is a parameter
        parfield = ['params_', tempparam];
        if isfield(par,parfield)
            str = par.(parfield);
            if regexp(str,'\Wstep\W')
                [tk1, str1] = gettoken(str, 'step');
                % string of the type 10ns step 6 ns
                tk2 = strtrim(gettoken(str1, ';'));
                [minval, unit] = kvgetvalue(tk1);
                step = kvgetvalue(tk2);
                arr = (0:asize-1)*step+minval;
            elseif regexp(str,'\Wlogto\W')
                [tk1, str1] = gettoken(str, 'logto');
                % string of the type 10ns logto 60 ns
                tk2 = strtrim(gettoken(str1, ';'));
                [minval, unit] = kvgetvalue(tk1);
                maxval = kvgetvalue(tk2);
                arr = logspace(log10(minval), log10(maxval),asize);
            elseif regexp(str,'\Wto\W')
                [tk1, str1] = gettoken(str, 'to');
                % string of the type 10ns to 60 ns
                tk2 = strtrim(gettoken(str1, ';'));
                [minval, unit] = kvgetvalue(tk1);
                maxval = kvgetvalue(tk2);
                arr = (0:1/(asize-1):1)*(maxval-minval)+minval;
            else
                % string of the type 10ns, 20ns, 30ns;
                [str1] = gettoken(str, ';');
                [tk1, str1] = gettoken(str1, ',');
                while ~isempty(tk1)
                    [arr(end+1),unit] = kvgetvalue(tk1);
                    [tk1, str1] = gettoken(str1, ',');
                    if isempty(tk1) && ~isempty(str1)
                        tk1 = str1; str1 = [];
                    end
                end
            end
        else
            str = par.(['aquisition_', tempparam]);
            arr = (0:asize-1)';
            unit = 's';
        end
        
        % Unit normalization
        switch unit
            case 'G',
            case 'K',
            case 's',
            otherwise
                umax = max(abs(arr));
                for kk = length(koeff):-1:1
                    if umax > koeff(kk)
                        uk = koeff(kk);
                        unit = [prefix(kk), unit];
                        arr = arr./uk;
                        break;
                    end
                end
        end
        
        res.(axislabel(counter)) = arr';
        res.([axislabel(counter), 'label']) = [sweepax{k}.var{1}, ', ',unit];
        counter = counter + 1;
    end
end
return

function [tk,rstr] = gettoken(istr, tok)

pos = strfind(istr, tok);

if isempty(pos)
    tk=istr;
    rstr='';
else
    tk=strtrim(istr(1:pos-1));
    rstr=strtrim(istr(pos+length(tok):end));
end
return

function res = safeget(strct, fld, deflt)
% function res = safeget(strct, fld, deflt)
% Returns the field of the structure or default
% value if field is absent.
% if default value is array then return value
% has not less elements than in this array
% (not in char case)

if isfield(strct, fld)
    res = strct.(fld);
    sd = size(deflt, 2);
    sr = size(res, 2);
    if sr < sd
        if iscell(deflt)
            [res{sr+1:sd}] = deal(deflt{sr+1:sd});
        elseif ~ischar(deflt)
            res(sr+1:sd) = deflt(sr+1:sd);
        end
        
    end
else
    res = deflt;
end
return

function [val, unit, pref, pref_val] = kvgetvalue(str)
% KVGETVALUE read string value in ci-standard units

% [val, str_unit, str_koefficient] = kvgetvalue(str)

prefix = ['p', 'n', 'u', 'm', 'k', 'M', 'G', 'T'];
koeff  = [1e-12, 1e-9, 1e-6, 1e-3, 1e3, 1e6, 1e9, 1e12];

% find substring that is a floating-point literal and convert to number
[validx1,validx2] = regexp(str,'[\-\+]?[0-9]*(\.[0-9]+)?');
val = str2double(str(validx1:validx2));

% extract substring that contains the unit
unit = strtrim(str(validx2+1:end));

pref = '';
pref_val = 1;
if length(unit) > 1 && ~isempty(unit)
  kk = strfind(prefix, unit(1));
  if ~isempty(kk)
    val = val * koeff(kk);
    unit = unit(2:end);
    pref = prefix(kk);
    pref_val = koeff(kk);
  end
end
return
