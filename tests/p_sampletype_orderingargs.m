function ok = test()

% Ordering functions with 1 or 3 input arguments are accepted, others are
% rejected

Exp.SampleFrame = [];
Exp.SampleRotation = [];
Exp.MolFrame = [];
Exp.CrystalSymmetry = '';
Opt = struct;

fcns = {@(beta) cos(beta).^2, @(alpha,beta,gamma) cos(beta).^2, ...
  @(varargin) 1, @(beta,gamma) cos(beta).^2};
accepted = [true true true false];

for k = numel(fcns):-1:1
  Exp.Ordering = fcns{k};
  try
    Exp_ = runprivate('p_sampletype',Exp,Opt);
    ok(k) = accepted(k) && nargin(Exp_.Ordering)~=1;
  catch err
    ok(k) = ~accepted(k) && contains(err.message,'Exp.Ordering');
  end
end

end
