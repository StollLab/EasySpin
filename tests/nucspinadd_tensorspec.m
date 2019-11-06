function err = test(opt)

% Test all combinations of tensor specifications (isotropic, axial, rhombic, full)
%-------------------------------------------------------------------------------

A1 = {[],1,[1 2],[1 2 3],rand(3)};
A2 = {6,[6 7],[6 7 8],rand(3)};

A1str = {'none','iso','axial','rhombic','full'};
A2str = {'iso','axial','rhombic','full'};

for i1 = 1:numel(A1)
  for i2 = 1:numel(A2)
    A1_ = A1{i1};
    A2_ = A2{i2};
    clear Sys0
    Sys0.S = 1/2;
    if ~isempty(A1{i1})
      Sys0.Nucs = '1H';
      Sys0.A = A1{i1};
      nSysNuc = 1;
    else
      nSysNuc = 0;
    end
    Sys2 = nucspinadd(Sys0,'13C',A2{i2});
    
    full = numel(A1_)==9 || numel(A2_)==9;
    if full
      ok(i1,i2) = all(size(Sys2.A)==[(nSysNuc+1)*3 3]);
    else
      d = max(numel(A1_),numel(A2_));
      if d==1
        ok(i1,i2) = numel(Sys2.A)==nSysNuc+1;
      else
        ok(i1,i2) = all(size(Sys2.A)==[nSysNuc+1,d]);
      end
    end
  end
end

err = any(~ok(:));
