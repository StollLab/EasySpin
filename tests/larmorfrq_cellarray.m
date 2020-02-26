function ok = test()

% cell array of isotopes

nu = larmorfrq({'1H','14N'},300);
ok = numel(nu)==2;
