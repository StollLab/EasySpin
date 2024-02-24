function ok = test()

% Test values of perpendicular dipolar couplings at 1 nm distance

r = [0;0;1];  % nm

Tee = diptensor(gfree,gfree,r);
Tee_perp = Tee(1,1);
ok(1) = areequal(Tee_perp,+52.0410,1e-4,'abs');

TeH = diptensor(gfree,'1H',r);
TeH_perp = TeH(1,1);
ok(2) = areequal(TeH_perp,-0.0790644,1e-6,'abs');
