function ok = test()

% Confirm that dipolar coupling along inter-spin vector has correct sign

r = [0;0;1];  % nm

% Two electrons: dipolar coupling should be negative
Tee = diptensor(gfree,gfree,r);
ok(1) = Tee(3,3)<0;

% Electron and proton: dipolar coupling should be positive
TeH = diptensor(gfree,'1H',r);
ok(2) = TeH(3,3)>0;
