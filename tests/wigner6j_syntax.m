function [err,data] = test(opt,olddata)

% Test 1: input syntax test
%======================================================
a0 = wigner6j([4,7/2,3/2,3/2,1,3]);      % 1 argument
a1 = wigner6j(4,7/2,3/2,3/2,1,3);        % 6 arguments
a2 = wigner6j([4,7/2,3/2],[3/2,1,3]);    % 2 arguments
a3 = wigner6j([4 3/2],[7/2 1],[3/2 3]);  % 3 arguments

err = ~areequal([a1 a2 a3],[a0 a0 a0],1e-12,'abs');
data = [];
