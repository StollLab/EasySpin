function ok = test()

% Strains on nuclear parameters are rewritten for each isotopologue

% Natural Cu: A strain scales with gn, eeqQ strain with Q, eta strain not
Sys.Nucs = 'Cu';
Sys.A = [50 50 500];
Sys.Q = [10 0.2];
Sys.StrainPars = {'A(3)','Q(1)','Q(2)'};
Sys.StrainFWHM = [30 1 0.1];
iso = isotopologues(Sys);
[~,gn] = nucdata({'63Cu','65Cu'});
[~,~,qm] = nucdata({'63Cu','65Cu'});
ok(1) = numel(iso)==2 && isequal(iso(1).StrainPars,{'A(1,3)','Q(1,1)','Q(1,2)'});
V1 = iso(1).StrainModes; V2 = iso(2).StrainModes;
ok(2) = areequal(V2(:,1),V1(:,1)*gn(2)/gn(1),1e-10,'rel');
ok(3) = areequal(V2(:,2),V1(:,2)*qm(2)/qm(1),1e-10,'rel');
ok(4) = areequal(V2(:,3),V1(:,3),1e-10,'rel');
ok(5) = ~isfield(iso,'StrainFWHM');

% Natural N with Q strain: no Q strain for 15N (I = 1/2)
Sys = struct('Nucs','N','A',[5 6 7],'Q',1,'StrainPars',{{'Q'}},'StrainFWHM',0.2);
iso = isotopologues(Sys);
Nucs = {iso.Nucs};
i14 = find(strcmp(Nucs,'14N'));
i15 = find(strcmp(Nucs,'15N'));
ok(6) = isequal(iso(i14).StrainPars,{'Q(1,1)'}) && isempty(iso(i15).StrainPars);
Exp = struct('mwFreq',9.5,'Range',[330 350]);
try
  pepper(Sys,Exp);
  ok(7) = true;
catch
  ok(7) = false;
end

% Strain on g only: fields are passed on unchanged
Sys = struct('Nucs','Cu','g',[2 2.2],'A',[50 500],'StrainPars',{{'g(2)'}},'StrainFWHM',0.01);
iso = isotopologues(Sys);
ok(8) = isequal(iso(1).StrainPars,{'g(2)'}) && isequal(iso(2).StrainFWHM,0.01);

% (row,col) reference into a 1 x nNuc row of isotropic A values
Sys = struct('Nucs','1H,Cu','A',[10 50],'StrainPars',{{'A(1,2)'}},'StrainFWHM',3);
iso = isotopologues(Sys);
ok(9) = isequal(iso(1).StrainPars,{'A(2,1)'});

% sigma and nn references: rows and pairs are remapped, pairs with a
% spin-0 isotope are dropped, nn widths scale with gn
Sys = struct('Nucs','1H,C,1H','A',[5 10 3],'nn',[0.1 0.2 0.3],'sigma',[1 1 1.001], ...
  'StrainPars',{{'nn(2)','sigma(3)','nn(3)'}},'StrainFWHM',[0.05 1e-4 0.05]);
iso = isotopologues(Sys);
i12 = strcmp({iso.Nucs},'1H,1H');
i13 = strcmp({iso.Nucs},'1H,13C,1H');
ok(10) = isequal(iso(i12).StrainPars,{'nn(1,1)','sigma(2,1)'});
ok(11) = isequal(iso(i13).StrainPars,{'nn(2,1)','sigma(3,1)','nn(3,1)'});
Sys = struct('Nucs','1H,Cu','A',[5 50],'nn',1,'StrainPars',{{'nn'}},'StrainFWHM',0.1);
iso = isotopologues(Sys);
ok(12) = areequal(iso(2).StrainModes,iso(1).StrainModes*gn(2)/gn(1),1e-12,'rel');
