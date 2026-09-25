function ok = test()

% Exp fields that conflict with Ori or B throw an error

Sys = struct('S',1,'D',[3000 500]);
conflictingFields = {'SampleFrame',[0 0 0]; 'Range',[0 600]};

hFig = figure('Visible','off');
for k = size(conflictingFields,1):-1:1
  Exp = struct('mwFreq',9.5,conflictingFields{k,1},conflictingFields{k,2});
  try
    levelsplot(Sys,'z',[0 600],Exp);
    ok(k) = false;
  catch err
    ok(k) = contains(err.message,conflictingFields{k,1});
  end
end
close(hFig);
