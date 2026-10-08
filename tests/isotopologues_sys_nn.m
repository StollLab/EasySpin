function ok = test()

% Isotopologues with nn and nnFrame: pairs with spin-0 isotopes are removed,
% and couplings are scaled with the gn ratios of both nuclei

% Spin-0 isotope in the middle: pairs (1,2),(1,3),(2,3) -> (1,3) only
Sys.Nucs = '1H,C,1H';
Sys.A = [5 10 3];
Sys.nn = [0.1 0.2 0.3];
Sys.nnFrame = [1 1 1; 2 2 2; 3 3 3];

Iso = isotopologues(Sys);
i12 = strcmp({Iso.Nucs},'1H,1H');
i13 = strcmp({Iso.Nucs},'1H,13C,1H');
ok(1) = isequal(Iso(i12).nn,0.2) && isequal(Iso(i12).nnFrame,[2 2 2]);
ok(2) = isequal(Iso(i13).nn,[0.1; 0.2; 0.3]) && isequal(Iso(i13).nnFrame,Sys.nnFrame);

% Scaling with gn, and zero coupling between nuclei of the same group
[~,gn] = nucdata({'63Cu','65Cu'});
r = gn(2)/gn(1);
Sys = struct('Nucs','1H,Cu','n',[1 2],'A',[5 50],'nn',1);
Iso = isotopologues(Sys);
i63 = strcmp({Iso.Nucs},'1H,63Cu');
i65 = strcmp({Iso.Nucs},'1H,65Cu');
iMix = strcmp({Iso.Nucs},'1H,63Cu,65Cu');
ok(3) = isequal(Iso(i63).nn,1);
ok(4) = areequal(Iso(i65).nn,r,1e-12,'rel');
ok(5) = areequal(Iso(iMix).nn,[1; r; 0],1e-12,'abs');

% Only one nucleus left: no couplings
Sys = struct('Nucs','1H,C','A',[5 10],'nn',0.1);
Iso = isotopologues(Sys);
i12 = strcmp({Iso.Nucs},'1H');
ok(6) = isempty(Iso(i12).nn);
