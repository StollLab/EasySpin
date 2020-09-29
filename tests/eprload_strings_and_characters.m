function ok = test()

FilePaths = {'./eprfiles/00012107','./eprfiles/00011201',"./eprfiles/00012107.dta","./eprfiles/00011201.spc"};

try
  for i = 1 : length(FilePaths)
    [~,~] = eprload(FilePaths{i});
  end
  ok = true;
catch
  ok = false;
end

