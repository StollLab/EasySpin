function ok = test()

% Check that invalid inputs raise errors

cases = {
  @() equivcouple(0.3,2),      'multiple of 1/2'
  @() equivcouple(-1/2,2),     'multiple of 1/2'
  @() equivcouple([1/2 1],2),  'multiple of 1/2'
  @() equivcouple('a',2),      'multiple of 1/2'
  @() equivcouple(1/2,2.5),    'nonnegative integer'
  @() equivcouple(1/2,-1),     'nonnegative integer'
  @() equivcouple(1/2,[2 3]),  'nonnegative integer'
  };

for k = size(cases,1):-1:1
  try
    cases{k,1}();
    ok(k) = false;
  catch err
    ok(k) = contains(err.message,cases{k,2});
  end
end

end
