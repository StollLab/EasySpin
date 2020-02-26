function ok = test()

% several fields and several nuclei

B = 100:150;
nu = larmorfrq('14N,15N,63Cu',100:150);
ok = all(size(nu)==[numel(B) 3]);
