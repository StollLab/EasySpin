function ok = test()

% Rotation matrices with beta close to 0 or pi that are rounded to a few
% decimals (e.g. tensor orientations printed by ORCA) have elements in
% R(3,1:2) and R(1:2,3) that are zero or pure rounding noise. alpha and gamma
% are then individually ill-defined, but the rotation matrix must still be
% recovered to within the rounding error.

% g-tensor orientations from ORCA output, columns are the principal axes
R_g2M{1} = [ 0.7704916  -0.6374502   0.0000000
             0.6374502   0.7704916  -0.0000000
            -0.0000000   0.0000001   1.0000000];
R_g2M{2} = [ 0.7704881   0.6374544  -0.0000000
             0.6374544  -0.7704881  -0.0000001
            -0.0000000   0.0000000  -1.0000000];
for k = 1:numel(R_g2M)
  R = R_g2M{k}.';
  ok(k) = norm(erot(eulang(R))-R) < 1e-6;
end

% Random matrices with beta close to 0 or pi, rounded to nDecimals
rng(4417);
nTrials = 200;
for t = 1:nTrials
  nDecimals = randi([7 12]);
  beta = pi*(rand<0.5) + 10^-randi([8 16])*randn;
  beta = min(abs(beta),2*pi-abs(beta));
  R = round(erot(2*pi*rand,beta,2*pi*rand)*10^nDecimals)/10^nDecimals;
  ok(end+1) = norm(erot(eulang(R))-R) < 10^(1-nDecimals); %#ok<AGROW>
end
