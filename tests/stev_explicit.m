function ok = test()

S = 9/2;

[Sz,Sp,Sm,I] = sop(S,'z','+','-','e');
s = S*(S+1);
cp = 1/2;
cm = 1/2i;

f = @(A,B) (A*B+B*A)/2;

% Rank 2
Op{1} = cp*(Sp^2+Sm^2);
Op{2} = cp*f(Sz,Sp+Sm);
Op{3} = 3*Sz^2 - s*I;
Op{4} = cm*f(Sz,Sp-Sm);
Op{5} = cm*(Sp^2-Sm^2);

Op0{1} = stev(S,[2,2]);
Op0{2} = stev(S,[2,1]);
Op0{3} = stev(S,[2,0]);
Op0{4} = stev(S,[2,-1]);
Op0{5} = stev(S,[2,-2]);

for k = 1:5
  ok2(k) = areequal(Op{k},Op0{k},1e-10,'abs');
end

% Rank 4
Op{1} = cp*(Sp^4+Sm^4);
Op{2} = cp*f(Sz,Sp^3+Sm^3);
Op{3} = cp*f(7*Sz^2-(s+5)*I,Sp^2+Sm^2);
Op{4} = cp*f(7*Sz^3-(3*s+1)*Sz,Sp+Sm);
Op{5} = 35*Sz^4-(30*s-25)*Sz^2+(3*s^2-6*s)*I;
Op{6} = cm*f(7*Sz^3-(3*s+1)*Sz,Sp-Sm);
Op{7} = cm*f(7*Sz^2-(s+5)*I,Sp^2-Sm^2);
Op{8} = cm*f(Sz,Sp^3-Sm^3);
Op{9} = cm*(Sp^4-Sm^4);

Op0{1} = stev(S,[4,4]);
Op0{2} = stev(S,[4,3]);
Op0{3} = stev(S,[4,2]);
Op0{4} = stev(S,[4,1]);
Op0{5} = stev(S,[4,0]);
Op0{6} = stev(S,[4,-1]);
Op0{7} = stev(S,[4,-2]);
Op0{8} = stev(S,[4,-3]);
Op0{9} = stev(S,[4,-4]);

for k = 1:9
  ok4(k) = areequal(Op{k},Op0{k},1e-10,'abs');
end

% Rank 6
Op{1} = cp*(Sp^6+Sm^6);
Op{2} = cp*f(Sz,Sp^5+Sm^5);
Op{3} = cp*f(11*Sz^2-(s+38)*I,Sp^4+Sm^4);
Op{4} = cp*f(11*Sz^3-(3*s+59)*Sz,Sp^3+Sm^3);
Op{5} = cp*f(33*Sz^4-(18*s+123)*Sz^2+(s^2+10*s+102)*I,Sp^2+Sm^2);
Op{6} = cp*f(33*Sz^5-(30*s-15)*Sz^3+(5*s^2-10*s+12)*Sz,Sp+Sm);
Op{7} = 231*Sz^6-(315*s-735)*Sz^4+(105*s^2-525*s+294)*Sz^2-(5*s^3-40*s^2+60*s)*I;
Op{8} = cm*f(33*Sz^5-(30*s-15)*Sz^3+(5*s^2-10*s+12)*Sz,Sp-Sm);
Op{9} = cm*f(33*Sz^4-(18*s+123)*Sz^2+(s^2+10*s+102)*I,Sp^2-Sm^2);
Op{10} = cm*f(11*Sz^3-(3*s+59)*Sz,Sp^3-Sm^3);
Op{11} = cm*f(11*Sz^2-(s+38)*I,Sp^4-Sm^4);
Op{12} = cm*f(Sz,Sp^5-Sm^5);
Op{13} = cm*(Sp^6-Sm^6);

Op0{1} = stev(S,[6,6]);
Op0{2} = stev(S,[6,5]);
Op0{3} = stev(S,[6,4]);
Op0{4} = stev(S,[6,3]);
Op0{5} = stev(S,[6,2]);
Op0{6} = stev(S,[6,1]);
Op0{7} = stev(S,[6,0]);
Op0{8} = stev(S,[6,-1]);
Op0{9} = stev(S,[6,-2]);
Op0{10} = stev(S,[6,-3]);
Op0{11} = stev(S,[6,-4]);
Op0{12} = stev(S,[6,-5]);
Op0{13} = stev(S,[6,-6]);

for k = 1:13
  ok6(k) = areequal(Op{k},Op0{k},1e-10,'abs');
end

ok = [ok2 ok4 ok6];
