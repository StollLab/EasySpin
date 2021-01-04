function ok = test()

% string('...') is used instead of "..." for compatibility with R2016b
FilePaths = {'./eprfiles/00012107','./eprfiles/00011201',...
  string('./eprfiles/00012107'),string('./eprfiles/00012107.dta'),...
  string('./eprfiles/00011201'),string('./eprfiles/00011201.spc')};

try
  for i = 1 : length(FilePaths)
    [~,~] = eprload(FilePaths{i});
  end
  ok = true;
catch
  ok = false;
end

