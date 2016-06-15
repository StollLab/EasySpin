function y = renorm(yi)
mi = 0;
ma = 1;

maxyi = max(yi);
minyi = min(yi);

y = mi + (ma-mi)*(yi-minyi)/(maxyi-minyi);

return
