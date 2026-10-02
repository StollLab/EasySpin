function ok = test()

% Check that invalid angles, options and output syntaxes raise errors

cases = {
  @() erot([1 2],0,0),         'finite real numbers'
  @() erot([],0,0),            'finite real numbers'
  @() erot(NaN,0,0),           'finite real numbers'
  @() erot(0,Inf,0),           'finite real numbers'
  @() erot(0,0,1i),            'finite real numbers'
  @() erot('abc'),             'finite real numbers'
  @() erot([1 2]),             'Three angles'
  @() erot([0 0 0],'diag'),    'either ''rows'' or ''cols'''
  @() erot([0 0 0],1),         'either ''rows'' or ''cols'''
  @() threeoutputs(),          'columns or the 3 rows'
  @() twooutputs(),            'Wrong number of outputs'
  @() oneoutputwithoption(),   '3 outputs required'
  };

for k = size(cases,1):-1:1
  try
    cases{k,1}();
    ok(k) = false;
  catch err
    ok(k) = contains(err.message,cases{k,2});
  end
end

% String option is accepted
[x,y,z] = erot([0.1 0.2 0.3],"rows");
ok(end+1) = areequal([x y z].',erot([0.1 0.2 0.3]),1e-12,'abs');

end

function threeoutputs()
[~,~,~] = erot(0,0,0);
end

function twooutputs()
[~,~] = erot(0,0,0);
end

function oneoutputwithoption()
R = erot(0,0,0,'cols'); %#ok<NASGU>
end
