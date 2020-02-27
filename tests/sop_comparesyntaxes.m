function ok = test()

% Compare the two different character-based syntaxes 
%
%  Syntax 1:  'eex'
%  Syntax 2:  'x3'

testfun = @(spins,syntax1,syntax2) areequal(sop(spins,syntax1),sop(spins,syntax2),1e-10,'abs');

ok(1) = testfun(1,'x','x1');
ok(2) = testfun([1/2 1],'xe','x1');
ok(3) = testfun([1/2 1],'xy','x1y2');
ok(4) = testfun([1/2 1],'e+','+2');
ok(5) = testfun([1/2 1 1/2],'ezz','z2z3');
ok(6) = testfun([1/2 1 1/2],'ezz','z3z2');
