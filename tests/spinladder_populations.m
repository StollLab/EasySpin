function ok = test()

J = 10e3;
Temp = 10;

S1 = 5/2;
S2 = 2;

Sys.S = [S1 S2];
Sys.ee = J;

C = spinladder(Sys,Temp);

for k=1:numel(C)
  weights(k) = C{k}.weight;
end

S  = abs(S1-S2):S1+S2;
E0 = J/2*(S.*(S+1)-S1*(S1+1)-S2*(S2+1));

% Compute energies and populations
pop = exp(-E0*1e6*planck/boltzm/Temp);
pop = pop.*(2*S+1);
pop = pop/sum(pop);

[E0,idx] = sort(E0);
pop = pop(idx);

ok = areequal(weights,pop,1e-10,'abs');
