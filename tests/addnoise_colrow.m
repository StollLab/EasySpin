function [err,data] = test(opt,olddata)

%======================================================
% Make sure addnoise works for column and row vectors without crashing
%======================================================

y = rand(100,1);
SNR = 10;

yf = addnoise(y,SNR,'f');
yn = addnoise(y,SNR,'f');
yu = addnoise(y,SNR,'f');

y = rand(1,100);
yf = addnoise(y,SNR,'f');
yn = addnoise(y,SNR,'f');
yu = addnoise(y,SNR,'f');

err = false;

data = [];
