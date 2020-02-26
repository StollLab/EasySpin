function ok = test()

% Empty isotope list

nu1 = larmorfrq('',150);
nu2 = larmorfrq('1H',[]);
ok = isempty(nu1) && isempty(nu2);
