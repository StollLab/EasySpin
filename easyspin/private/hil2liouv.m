% hil2liouv takes the rank-0, rank-1, and rank-2 rotational basis operators
% in Hilbert space and transforms them into Liouville space

function [Q0L,Q1L,Q2L] = hil2liouv(Q0H,Q1H,Q2H)

nI = length(Q0H);
kronkron = @(A) speyekron(nI,A)-spkroneye(A.',nI);

% rank-0
Q0L = kronkron(Q0H);

% rank-1
if ~isempty(Q1H)
  Q1L = cell(3,3);
  for im1 = 1:3
    for im2 = 1:3
      Q1L{im1,im2} = kronkron(Q1H{im1,im2});
    end
  end
else
  Q1L = {};
end

% rank-2
Q2L = cell(5,5);
for im1 = 1:5
  for im2 = 1:5
    Q2L{im1,im2} = kronkron(Q2H{im1,im2});
  end
end

return