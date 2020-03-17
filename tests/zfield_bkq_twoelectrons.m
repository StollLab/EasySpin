function ok = test()

% Test Bkq terms

B4a = rand(1,9);
B4b = rand(1,9);
clear Sys
Sys.S = [7/2 7/2];
Sys.B4 = [B4a; B4b];
Sys.ee = 1e-70;

% (1) Calculate rank-4 interaction using zfield
Op1 = zfield(Sys);

% (2) Calculate rank-4 interaction manually
Op2 = 0;
k = 4;
q = k:-1:-k;
for iq = 1:numel(q)
  for iSpin = 1:2
    Op2 = Op2 + Sys.B4(iSpin,iq)*stev(Sys.S,[k,q(iq),iSpin]);
  end
end

ok = areequal(Op1,Op2,1e-10,'abs');
