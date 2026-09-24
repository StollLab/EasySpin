function ok = test()

% Check that the hyperfine cutoff removes exactly the nuclei with hyperfine
% couplings at or below the cutoff, and that all per-nucleus fields stay
% consistent.

files = {
  'v4.2.1/nitroxide.prop'
  'v6.1.0/nitroxide.out'
  'v6.1.0/nitroxide.property.txt'
  };
cutoffs = [0 5 20];  % MHz

perNucFields = {'A','AFrame','Q','QFrame'};

ok = true(numel(files),numel(cutoffs));
for f = 1:numel(files)
  Sys0 = orca2easyspin(['orca/' files{f}]);
  Nucs0 = strsplit(Sys0.Nucs,',');
  for c = 1:numel(cutoffs)
    Sys = orca2easyspin(['orca/' files{f}],cutoffs(c));
    keep = max(abs(Sys0.A),[],2).' > cutoffs(c);

    okk = ischar(Sys.Nucs) && strcmp(Sys.Nucs,strjoin(Nucs0(keep),','));
    okk = okk && isequal(Sys.NucsIdx,Sys0.NucsIdx(keep));
    for p = perNucFields
      okk = okk && isequal(Sys.(p{1}),Sys0.(p{1})(keep,:));
    end
    ok(f,c) = okk && any(~keep)==(cutoffs(c)>0);
  end
end

ok = ok(:);
