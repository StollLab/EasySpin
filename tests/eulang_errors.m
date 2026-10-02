function ok = test()

% Check that invalid input matrices and output syntaxes raise errors

R_approx = [ 0.0   1.0  -0.1
            -0.8   0.0   0.6
             0.6   0.1   0.8];

cases = {
  @() eulang(eye(2)),              'real-valued 3x3 matrix'
  @() eulang(eye(3,4)),            'real-valued 3x3 matrix'
  @() eulang(ones(3,3,2)),         'real-valued 3x3 matrix'
  @() eulang(1i*eye(3)),           'real-valued 3x3 matrix'
  @() eulang(nan(3)),              'real-valued 3x3 matrix'
  @() eulang(diag([Inf 1 1])),     'real-valued 3x3 matrix'
  @() eulang(R_approx),            'not orthogonal'
  @() eulang(-eye(3)),             'negative determinant'
  @() eulang(diag([1 1 -1])),      'negative determinant'
  @() twooutputs(eye(3)),          'Wrong number of output arguments'
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

function twooutputs(R)
[~,~] = eulang(R);
end
