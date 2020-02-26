function ok = test()

%======================================================
% Make sure addnoise works for column and row vectors without crashing
%======================================================

rng(454);

ycol = rand(100,1);
SNR = 10;

yf = addnoise(ycol,SNR,'f');
yn = addnoise(ycol,SNR,'n');
yu = addnoise(ycol,SNR,'u');

yrow = ycol.';
yf = addnoise(yrow,SNR,'f');
yn = addnoise(yrow,SNR,'n');
yu = addnoise(yrow,SNR,'u');

ok = true;
