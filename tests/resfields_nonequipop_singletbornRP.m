function ok = test()

% Test correct sign of emission/absorption for spin-correlated radical pair
% EAEA polarization for a singlet-born radical pair with J<0 in field-swept spectrum

Sys.S = [1/2 1/2];
Sys.g = [2.00 1.99];
Sys.J = -5; % MHz

S = 1/sqrt(2)*[0; 1; -1; 0];
Sys.initState = S*S';

% Field sweep
Exp = struct('mwFreq',9.75,'Range',[345 353]);
Exp.Harmonic = 0;

[B,A1] = resfields(Sys,Exp);
[B,idx] = sort(B);
A1 = A1(idx);

% Frequency sweep
Expf = struct('Field',349);

[f,A2] = resfreqs_matrix(Sys,Expf);
[f,idx] = sort(f);
A2 = A2(idx);

if numel(B)==4 && numel(f)==4
  ok = (A1(1)<0 && A1(2)>0 && A1(3)<0 && A1(4)>0) && ...
       (A2(1)>0 && A2(2)<0 && A2(3)>0 && A2(4)<0);
else
  ok = false;
end
