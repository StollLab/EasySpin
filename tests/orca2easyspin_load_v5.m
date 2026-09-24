function ok = test()

% Load all supported ORCA 5 file types for all test molecules, and check
% the structure and contents of the spin system.

files = {
  'v5.0.4/%s.out'
  };

%          name                   S    elements           has D
mols = {
  'aminoxyl',             1/2, 'N,O,H,H',              false
  'dioxygen',             1,   'O,O',                  false
  'nitroxide',            1/2, 'C,C,C,N,C,O,H,H,H,H,H,H', false
  'tripletformaldehyde',  1,   'H,C,H,O',              true
  };

ok = true(numel(files),size(mols,1));
for iFile = 1:numel(files)
  for iMol = 1:size(mols,1)
    [name,S,elements,hasD] = mols{iMol,:};
    Sys = orca2easyspin(['orca/' sprintf(files{iFile},name)]);

    Elements = strsplit(elements,',');
    nAtoms = numel(Elements);
    nNucs = numel(Sys.NucsIdx);

    okk = isscalar(Sys) && Sys.S==S && Sys.charge==0;
    okk = okk && isequal(Sys.Elements,Elements);
    okk = okk && isequal(size(Sys.xyz),[nAtoms 3]);
    okk = okk && isequal(size(Sys.g),[1 3]) && isequal(size(Sys.gFrame),[1 3]);
    okk = okk && isfield(Sys,'D')==hasD && isfield(Sys,'DFrame')==hasD;
    if hasD
      okk = okk && isequal(size(Sys.D),[1 3]) && isequal(size(Sys.DFrame),[1 3]);
    end

    % All atoms are computed, all have non-zero hyperfine coupling
    okk = okk && ischar(Sys.Nucs) && strcmp(Sys.Nucs,elements);
    okk = okk && isequal(Sys.NucsIdx,1:nAtoms);
    for f = {'A','AFrame','Q','QFrame'}
      okk = okk && isequal(size(Sys.(f{1})),[nNucs 3]);
    end
    okk = okk && all(any(Sys.A,2));

    ok(iFile,iMol) = okk;
  end
end

ok = ok(:);
