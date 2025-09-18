function ok = test()

gyromagratios = gamman('1H,2H,15N');

ok = numel(gyromagratios)==3;

end
