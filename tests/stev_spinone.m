function [err,data] = test(opt,olddata)

% Test 2: all operators for S=1
%======================================================
a = sqrt(2)/2;
Op{1} = 1i*[0 -a 0; a 0 -a; 0 a 0];
Op{2} = [1 0 0; 0 0 0; 0 0 -1];
Op{3} = [0 a 0; a 0 a; 0 a 0];

clear ok
q = -1:1;
for iq=1:numel(q)
  Opp = stev(1,1,q(iq));
  ok(iq) = all(abs(Opp(:)-Op{iq}(:))<1e-10);
end

err = any(~ok);
data = [];
