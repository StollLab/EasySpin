function [err,data] = test(opt,olddata)
%==================================================================
% Test whether the tensor output of zeemanho is identical to the full
% matrix output
%==================================================================


Sys.S = randi(10)/2+4;

Sys.ZB02.l =2;
Sys.ZB04.l =4;
Sys.ZB06.l =6;
Sys.ZB08.l =8;

Sys.ZB11.l= [0,2];

Sys.ZB22.l = 0:2:4;
Sys.ZB24.l = 2:2:6;
Sys.ZB26.l = 4:2:8;
Sys.ZB28.l = 6:2:10;

Sys.ZB31.l= 2:2:4;
Sys.ZB33.l= 0:2:6;
Sys.ZB35.l= 2:2:8;
Sys.ZB37.l= 4:2:10;

%fill list of coefficients with random numbers
strlist =fieldnames(Sys);
for n = length(strlist):-1:1
  if ~ strncmp(strlist{n}, 'ZB',2), continue; end;
  for m = length(Sys.(strlist{n}).l):-1:1  
    Sys.(strlist{n}).vals{m} = rand(2*Sys.(strlist{n}).l(m)+1,1);
  end
end

%obtain tensors
[G0,G1,G2,G3] = zeemanho(Sys);
zHo = zeemanho(Sys);

%obtain full matrices for the corresponding order
B = rand(3,1);
H{1} = zeemanho(Sys,B, [],'',0);
H{2} = zeemanho(Sys,B, [],'',1);
H{3} = zeemanho(Sys,B, [],'',2);
H{4} = zeemanho(Sys,B, [],'',3);

%build full matrices from tensors
size = 2*Sys.S+1;
for n = 4:-1:1
  Ht1{n} = zeros(size);
  Ht2{n} = zeros(size);
end
Ht1{1} = G0;
Ht2{1} = zHo{1};
for n=1:3
  Ht1{2} = Ht1{2}+ G1{n}*B(n);
  Ht2{2} = Ht2{2}+ zHo{2}{n}*B(n);  
  for m=1:3
    Ht1{3} = Ht1{3}+ G2{n,m}*B(n)*B(m);
    Ht2{3} = Ht2{3}+ zHo{3}{n,m}*B(n)*B(m);
    for o=1:3
      Ht1{4} = Ht1{4}+ G3{n,m,o}*B(n)*B(m)*B(o);
      Ht2{4} = Ht2{4}+ zHo{4}{n,m,o}*B(n)*B(m)*B(o);      
    end
  end
end

% test
threshold = 1e-4;
for n = 4:-1:1
 errl(n) = ~areequal(H{n},Ht1{n},threshold);
 errl(4+n) = ~areequal(H{n},Ht2{n},threshold);
end
err = any(errl);
data =[];
