function [err,data] = test(opt,olddata)
%==================================================================
% Test whether the tensor output of zeemanho is identical to the full
% matrix output
%==================================================================

rng_(5,'twister'); % seed randon number generator

Sys.S = 9/2;

for lB = 0:3
  for lS = 1:8
    if mod((lB+lS),2) == 0
      for l = abs(lB-lS): (lB+lS)
        str = ['Ham',num2str([lB,lS,l],'%i%i%i')];
        Sys.(str) = [];
      end
    end
  end
end

% Fill list of coefficients with random numbers
strlist = fieldnames(Sys);
for n = length(strlist):-1:1
  if ~strncmp(strlist{n},'Ham',3), continue; end
  l = str2num(strlist{n}(6:end));
  Sys.(strlist{n}) = rand(1,2*l+1); 
end

% Obtain tensors
[G0,G1,G2,G3] = zeemanho(Sys);
zHo = zeemanho(Sys);

% Obtain full matrices for the corresponding order
B = rand(3,1);
H{1} = zeemanho(Sys,B, [],'sparse',0);
H{2} = zeemanho(Sys,B, [],'sparse',1);
H{3} = zeemanho(Sys,B, [],'sparse',2);
H{4} = zeemanho(Sys,B, [],'sparse',3);

% Build full matrices from tensors
size = 2*Sys.S+1;
for n = 4:-1:1
  Ht1{n} = zeros(size);
  Ht2{n} = zeros(size);
end
Ht1{1} = G0;
Ht2{1} = zHo{1};
for n = 1:3
  Ht1{2} = Ht1{2}+ G1{n}*B(n);
  Ht2{2} = Ht2{2}+ zHo{2}{n}*B(n);  
  for m = 1:3
    Ht1{3} = Ht1{3}+ G2{n,m}*B(n)*B(m);
    Ht2{3} = Ht2{3}+ zHo{3}{n,m}*B(n)*B(m);
    for o = 1:3
      Ht1{4} = Ht1{4}+ G3{n,m,o}*B(n)*B(m)*B(o);
      Ht2{4} = Ht2{4}+ zHo{4}{n,m,o}*B(n)*B(m)*B(o);      
    end
  end
end

% test
threshold = 1e-4;
err = false;
for n = 1:4
 err = err || ~areequal(H{n},Ht1{n},threshold,'abs') || ~areequal(H{n},Ht2{n},threshold,'abs');
end
data =[];
