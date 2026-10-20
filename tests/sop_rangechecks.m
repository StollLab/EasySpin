function ok = test()

% Check that out-of-range spin indices and levels throw errors

badInputs = {
  {[1/2 1/2],'x1,z5'}      % spin index too large, comma syntax
  {[1/2 1/2],'x(1|2)3'}    % spin index too large, level syntax
  {[1/2 1/2],'x0,z1'}      % spin index zero
  {[1/2 1/2],[1 1; 5 3]}   % spin index too large, numeric syntax
  {[1/2 1/2],[0 1]}        % spin index zero, numeric syntax
  {1,'x(1|4)'}             % level too large
  {1,'x(0|1)'}             % level zero
  {[1/2 1],'e(3)1,z2'}     % level too large for spin 1
  };

for k = 1:numel(badInputs)
  try
    sop(badInputs{k}{:});
    ok(k) = false;
  catch
    ok(k) = true;
  end
end

% Valid edge cases still work
ok(end+1) = isequal(size(sop([1/2 1],'e(2)1,x(2|3)2')),[6 6]);
ok(end+1) = isequal(size(sop([1/2 1],[2 1])),[6 6]);
