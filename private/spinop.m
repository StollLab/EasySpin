% spinop takes an array of spins as input and forms the x, y, and z components of
% the spin vector operators for each spin in the array. The vector operator
% components are stored in a cell array in which each row corresponds to a
% spin and the three columns are the x, y, and z components respectively.
% The output operators are in the direct product space of all of the spins.

function SpinOps = spinop(Spin)

nSpins = numel(Spin);


for iSpin = 1:nSpins
  SpinOps{iSpin,1} = sop(Spin,iSpin,1);
  SpinOps{iSpin,2} = sop(Spin,iSpin,2);
  SpinOps{iSpin,3} = sop(Spin,iSpin,3);
end

 
return