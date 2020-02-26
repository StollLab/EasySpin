function ok = test()

% explicit comparison

Sys = struct('S',[1/2],'Nucs','1H','g',[2 2 3]);
Sys.A = rand(1,3);
Sys.AFrame = rand(1,3)*2*pi;

R = erot(Sys.AFrame);
A = R'*diag(Sys.A)*R;
HFI = 0;
for k=1:3
  for q=1:3
    HFI = HFI + A(k,q)*sop(Sys,[1 2],[k q]);
  end
end
H1 = hfine(Sys);

ok = areequal(HFI,H1,1e-10,'rel');
