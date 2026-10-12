function ok = test()

% Input forms of Sys.lwAnisoG and its older name Sys.HStrain

Sys0 = struct('S',1/2,'g',[2 2.1 2.2]);

warnState = warning('off','EasySpin:HStrainRenamed');

% accepted inputs and their normalized forms
Inputs = {[], 5, [5 7], [5 6 7]};
Results = {[0 0 0], [5 5 5], [5 5 7], [5 6 7]};

% field absent
Sys = runprivate('validatespinsys',Sys0);
ok = isequal(Sys.lwAnisoG,[0 0 0]);

Names = {'lwAnisoG','HStrain'};
for n = 1:numel(Names)
  for k = 1:numel(Inputs)
    Sys_ = Sys0;
    Sys_.(Names{n}) = Inputs{k};
    lastwarn('');
    [Sys,err] = runprivate('validatespinsys',Sys_);
    ok(end+1) = isempty(err) && isequal(Sys.lwAnisoG,Results{k}) && ~isfield(Sys,'HStrain'); %#ok<AGROW>
    % warning only for the older name
    [~,warnId] = lastwarn;
    ok(end+1) = strcmp(warnId,'EasySpin:HStrainRenamed')==strcmp(Names{n},'HStrain'); %#ok<AGROW>
  end

  % rejected inputs: column vector, too many elements, negative value
  BadInputs = {[5;6;7], [5 6 7 8], [5 -6 7]};
  for k = 1:numel(BadInputs)
    Sys_ = Sys0;
    Sys_.(Names{n}) = BadInputs{k};
    [~,err] = runprivate('validatespinsys',Sys_);
    ok(end+1) = contains(err,'lwAnisoG'); %#ok<AGROW>
  end
end

% both names given
Sys_ = Sys0;
Sys_.lwAnisoG = [5 6 7];
Sys_.HStrain = [5 6 7];
[~,err] = runprivate('validatespinsys',Sys_);
ok(end+1) = contains(err,'HStrain');

warning(warnState);
