function newsop = expandsops(System,spinop,position)

Sys = System;
newsop = spinop;

if position == 1
  
  for iOtherSpin = 2:numel(Sys.Spins)
    I = eye(2*Sys.Spins(iOtherSpin)+1);
    newsop = kron(newsop,I);
  end
  
end


return