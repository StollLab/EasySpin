function ok = test()

% Test whether eulang works with all supported output syntaxes

angles = deg2rad([68 39 213]);
R = erot(angles);

eulang(R);
[~,~,~] = eulang(R);
[~] = eulang(R);

ok = true;
