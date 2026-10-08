function ok = test()

% Parsing of Sys.StrainPars: accepted references and errors

Sys.S = 1;
Sys.g = [2 2.1 2.2];
Sys.gFrame = [0.1 0.2 0.3];
Sys.D = [300 50];
Sys.Nucs = '14N,1H';
Sys.A = [10 12 14; 3 4 5];
Sys.Q = [2 0.3; 0 0];

% Accepted references, and all input forms give the same result
forms = {{'g(3)','D(1,2)','A(2,1)','Q(1,2)','gFrame(2)'}, ...
         {"g(3)","D(1,2)","A(2,1)","Q(1,2)","gFrame(2)"}, ...
         {'g(3)',"D(1,2)",'A(2,1)',"Q(1,2)",'gFrame(2)'}, ...
         ["g(3)" "D(1,2)" "A(2,1)" "Q(1,2)" "gFrame(2)"]};
Sys.StrainFWHM = [0.01 10 1 0.1 0.05];
for k = 1:numel(forms)
  Sys.StrainPars = forms{k};
  [Sys_,err] = runprivate('validatespinsys',Sys);
  ok(k) = isempty(err) && numel(Sys_.StrainData.Deriv)==5;
  if k==1, Deriv1 = Sys_.StrainData.Deriv; else, ok(k) = ok(k) && isequal(Deriv1,Sys_.StrainData.Deriv); end
end
Sys.StrainFWHM = 0.01;
Sys.StrainPars = 'g(3)';
[~,err] = runprivate('validatespinsys',Sys); ok(end+1) = isempty(err);
Sys.StrainPars = "g(1,3)";
[~,err] = runprivate('validatespinsys',Sys); ok(end+1) = isempty(err);

% Errors
bad = {'gStrain(1)','not supported';  % disallowed field (also obsolete name)
       'B2(1)','not supported';
       'g_(1)','not supported';
       'ee(1)','not given';             % missing field
       'g(4)','out of range';
       'g(2,1)','out of range';
       'g','more than one element';
       'g(1','not a valid reference';
       'Q(2,1)','spin 0'};              % nucleus with I=0? (1H: I=1/2 -> Q has no effect)
for k = 1:size(bad,1)
  Sys.StrainPars = bad(k,1);
  [~,err] = runprivate('validatespinsys',Sys);
  ok(end+1) = ~isempty(err); %#ok<*AGROW>
end

% Duplicates after resolving the index
Sys.StrainPars = {'g(3)','g(1,3)'}; Sys.StrainFWHM = [0.01 0.01];
[~,err] = runprivate('validatespinsys',Sys);
ok(end+1) = contains(err,'same parameter');

% Full 3x3 tensors: elements and frames rejected; 1x6 elements accepted
Sys2.g = [2 0 0; 0 2.1 0; 0 0 2.2];
Sys2.gFrame = [0 0 0];
Sys2.StrainFWHM = 0.01;
Sys2.StrainPars = {'g(1,1)'};
[~,err] = runprivate('validatespinsys',Sys2); ok(end+1) = contains(err,'full');
Sys2.StrainPars = {'gFrame(1)'};
[~,err] = runprivate('validatespinsys',Sys2); ok(end+1) = contains(err,'full');
Sys2.g = [2 2.1 2.2 0.01 0 0];
[~,err] = runprivate('validatespinsys',Sys2); ok(end+1) = contains(err,'full');
Sys2.StrainPars = {'g(4)'};
[~,err] = runprivate('validatespinsys',Sys2); ok(end+1) = isempty(err);

% Parameters without effect
Sys3 = struct('S',1/2,'D',100,'StrainPars',{{'D'}},'StrainFWHM',10);
[~,err] = runprivate('validatespinsys',Sys3); ok(end+1) = contains(err,'no effect');
Sys3 = struct('Nucs','1H','A',[5 6],'Q',0,'StrainPars',{{'Q'}},'StrainFWHM',1);
[~,err] = runprivate('validatespinsys',Sys3); ok(end+1) = contains(err,'no effect');
Sys3 = struct('Nucs','12C,1H','A',[5; 6],'StrainPars',{{'A(1)'}},'StrainFWHM',1);
[~,err] = runprivate('validatespinsys',Sys3); ok(end+1) = contains(err,'spin 0');
Sys3 = struct('g',2,'gFrame',[0.1 0.2 0.3],'StrainPars',{{'gFrame(2)'}},'StrainFWHM',0.1);
[~,err] = runprivate('validatespinsys',Sys3); ok(end+1) = contains(err,'no effect');

% ee strains need ee input, J strains need J input
Sys4 = struct('S',[1/2 1/2],'g',[2 2],'J',100,'StrainPars',{{'ee(1)'}},'StrainFWHM',1);
[~,err] = runprivate('validatespinsys',Sys4); ok(end+1) = contains(err,'not given');
Sys4 = struct('S',[1/2 1/2],'g',[2 2],'ee',100,'StrainPars',{{'J(1)'}},'StrainFWHM',1);
[~,err] = runprivate('validatespinsys',Sys4); ok(end+1) = contains(err,'not given');

% Malformed field: regular size error from validatespinsys
Sys5 = struct('Nucs','1H','A',[1 2 3 4],'StrainPars',{{'A(1)'}},'StrainFWHM',1);
[~,err] = runprivate('validatespinsys',Sys5); ok(end+1) = contains(err,'Size of Sys.A');

% Strain fields without StrainPars
Sys6 = struct('g',2,'StrainFWHM',0.01);
[~,err] = runprivate('validatespinsys',Sys6); ok(end+1) = contains(err,'StrainPars');
Sys6 = struct('g',2,'StrainPars',{{}},'StrainCorr',[]);
[~,err] = runprivate('validatespinsys',Sys6); ok(end+1) = contains(err,'StrainPars');
