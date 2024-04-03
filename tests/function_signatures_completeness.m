function ok = test(opt)

% Read and parse function signatures JSON file
esPath = fileparts(which('sop'));
fname =  [esPath filesep 'functionSignatures.json'];
fid = fopen(fname);
raw = fread(fid,Inf);
fclose(fid);
str = char(raw');
val = jsondecode(str);
funcsigs = fieldnames(val);

% Remove duplicates that end in _1, _2, _3, etc.
for i = numel(funcsigs):-1:1
  tail = funcsigs{i}(end-1:end);
  num = str2double(tail(2));
  if tail(1)=='_' && ~isempty(num), funcsigs(i) = []; end
end

% Get list of function files
funcs = dir([esPath filesep '*.m']);
funcs = {funcs.name}.';

% Remove .m from filename
for i = numel(funcs):-1:1
  funcs{i} = funcs{i}(1:end-2);
end

% Make joint list of functions and signatures, exclude functions where
% no signature information is provided
items = union(funcs,funcsigs);
exclusions = {'Contents', 'runprivate', 'eschecker',...
  'crystalfield', 'eeint', 'hfine', 'nnint', 'nquad', 'sham', ...
  'soint', 'zeeman', 'zeemanho', 'zfield','eigfields',...
  'esfit_genetic','esfit_grid','esfit_legacy','esfit_levmar','esfit_montecarlo',...
  'esfit_simplex','esfit_swarm','easyspinversioncheck','easyspinupdate'};  % items to exclude
items = setdiff(items,exclusions);

% Determine for each item whether function and signature are present
for i = numel(items):-1:1
  sig_present(i) = ~isempty(find(matches(funcsigs,items{i}),1));
  func_present(i) = ~isempty(find(matches(funcs,items{i}),1));
end

ok = all(sig_present & func_present);

if opt.Display
  missingFunctions = items(sig_present & ~func_present);
  missingSignatures = items(~sig_present & func_present);
  missingFunctions
  missingSignatures
  fprintf('%d missing functions, %d missing signatures\n',numel(missingFunctions),numel(missingSignatures));
end

end
