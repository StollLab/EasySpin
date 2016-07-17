function [err,data] = test(opt,olddata)

% Expansion of a single angle to a vector

N = 100;
n = ang2vec(rand,rand(1,N));
n = ang2vec(rand(1,N),rand);
n = ang2vec(rand,rand(N,1));
n = ang2vec(rand(N,1),rand);

err = 0;

data = [];
