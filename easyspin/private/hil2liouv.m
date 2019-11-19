% hil2liouv takes Hilbert-space operators and transforms them to Liouville space

function [HL] = hil2liouv(H)

kronkron = @(A) eyekron(length(A),A) - kroneye(A.',length(A));

if ~iscell(H)
  
  if ~issparse(H)
    error('Input matrix must be in sparse format.');
  end
  
  HL = kronkron(H);
  
else
  
  HL = cell(size(H));
  for k = 1:numel(H)
    if ~issparse(H{k})
      error('Input matrix must be in sparse format.');
    end
    HL{k} = kronkron(H{k});
  end
  
end
