function ok = test()

% Sys.lwAnisoG and its older name Sys.HStrain give the same orientation selection

Exp = struct('mwFreq',9.5,'Field',320,'ExciteWidth',100);
Opt = struct('GridSize',11);

Inputs = {50, [40 50 60]};
for k = 1:numel(Inputs)
  Sys = struct('S',1/2,'g',[2 2.1 2.2],'lwAnisoG',Inputs{k});
  w1 = orisel(Sys,Exp,Opt);
  Sys = struct('S',1/2,'g',[2 2.1 2.2],'HStrain',Inputs{k});
  warnState = warning('off','EasySpin:HStrainRenamed');
  w2 = orisel(Sys,Exp,Opt);
  warning(warnState);
  ok(k) = isequal(w1,w2);
end

% scalar is the same as three equal values
Sys = struct('S',1/2,'g',[2 2.1 2.2],'lwAnisoG',[50 50 50]);
w3 = orisel(Sys,Exp,Opt);
Sys = struct('S',1/2,'g',[2 2.1 2.2],'lwAnisoG',50);
w4 = orisel(Sys,Exp,Opt);
ok(end+1) = isequal(w3,w4);
