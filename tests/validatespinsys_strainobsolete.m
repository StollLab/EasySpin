function ok = test()

% Obsolete strain fields give an error with working replacement code

c = 0;
c=c+1; S{c} = struct('g',[2 2.1 2.2],'gStrain',[0.01 0 0.03]);
       expected{c} = {'g(1,1)','g(1,3)'};
c=c+1; S{c} = struct('g',[2 2.2],'gStrain',[0.01 0.02]);
       expected{c} = {'Sys.g = [2 2 2.2]','g(1,3)'};
c=c+1; S{c} = struct('Nucs','63Cu','g',[2 2.2],'A',[50 500],'gStrain',[0 0.02],'AStrain',[0 30]);
       expected{c} = {'A(1,3)','Sys.StrainCorr'};
c=c+1; S{c} = struct('Nucs','63Cu','g',[2 2.2],'A',[50 500],'gStrain',[0 0.02],'AStrain',[0 30],'gAStrainCorr',-1);
       expected{c} = {'-1'};
c=c+1; S{c} = struct('S',1,'D',[300 50],'DStrain',[30 10],'DStrainCorr',0.5);
       expected{c} = {'D(1,1)','D(1,2)','0.5'};
c=c+1; S{c} = struct('S',1,'D',300,'DStrain',[30 10]);
       expected{c} = {'Sys.D = [300 0]','D(1,2)'};
c=c+1; S{c} = struct('S',1,'D',300,'DStrain',30);
       expected{c} = {'D(1,1)'};
c=c+1; S{c} = struct('S',1/2,'gStrain',[]);
       expected{c} = {'no longer supported'};

for k = 1:c
  [~,err] = runprivate('validatespinsys',S{k});
  ok(k) = contains(err,'no longer supported');
  for e = 1:numel(expected{k})
    ok(k) = ok(k) && contains(err,expected{k}{e});
  end
  % Replacement code must give a valid spin system
  lines = regexp(err,'Sys\.[^\n]*','match');
  if isempty(lines), continue; end
  Sys = rmfield(S{k},intersect(fieldnames(S{k}),{'gStrain','AStrain','gAStrainCorr','DStrain','DStrainCorr'}));
  for L = 1:numel(lines)
    if startsWith(lines{L},'Sys.gStrain') || contains(lines{L},'no longer'), continue; end
    eval(lines{L});
  end
  [~,err2] = runprivate('validatespinsys',Sys);
  ok(k) = ok(k) && isempty(err2);
end
