% orca2easyspin_maintxt   Read EPR properties from main ORCA output file

function Sys = orca2easyspin_proptxt(propfilename,HyperfineCutoff)

% Read entire file into cell array
fh = fopen(propfilename);
allLines = textscan(fh,'%s','whitespace','','delimiter','\n');
allLines = allLines{1};
fclose(fh);

if ~contains(allLines{2},"!PROPERTIES!")
  error('This is not a valid text-based ORCA property file.');
end

% Locate lines with section titles
sections = find(cellfun(@(L)L(1)=='$',allLines));

S = [];
g = [];
gFrame = [];
D = [];
DFrame = [];
A = [];
AFrame = [];
Q = [];
QFrame = [];
E = [];  % total energy, hartree
charge = [];  % total charge, e
atoms = [];  % list of atoms (element symbols)
xyz = [];  % atom coordinates (in ångstrom)

% Loop over all sections and parse relevant information
for iSection = 1:numel(sections)
  idx0 = sections(iSection);
  sectionTitle = allLines{idx0};
  switch sectionTitle
    case '$ SCF_Energy'
      E = readvalue(allLines{idx0+4});

      if allLines{idx0+6}(1)~='$'
        error('Section SCF_Energy is longer than expected. Abording.');
      end

    case '$ DFT_Energy'
      % not relevant

    case '$ Calculation_Info'
      multiplicity = readvalue(allLines{idx0+4});
      S = (multiplicity-1)/2;
      charge = readvalue(allLines{idx0+5});

    case '$ SCF_Electric_Properties'
      % not relevant

    case '$ EPRNMR_GTensor'
      g_raw = readmatrix(allLines(idx0+(8:10)));
      g_vecs = readmatrix(allLines(idx0+(13:15)));
      g_vals = readvec(allLines{idx0+18}(9:end));
      g = g_vals;
      gFrame = eulang(g_vecs);

    case '$ EPRNMR_DTensor'
      % Bug in ORCA 5.0.2: D_vals is all-zero, even though D_raw is correct
      D_raw = readmatrix(allLines(idx0+(8:10)));
      %D_vals = readvec(allLines{idx0+13}(9:end));
      %D_vecs = readmatrix(allLines(idx0+(16:18)));
      [D_vecs,D_vals] = eig(D_raw);
      D_vals = diag(D_vals).';
      D = D_vals;
      DFrame = eulang(D_vecs);

    case '$ EPRNMR_ATensor'
      % colon is missing after "Number of stored nuclei"
      valuestr = regexp(allLines{idx0+4},'\d+$','match','once');
      nNuclei = str2double(valuestr);
      A = zeros(nNuclei,3);
      AFrame = zeros(nNuclei,3);
      i = idx0 + 7;
      for n = 1:nNuclei
        A_raw = readmatrix(allLines(i+(6:8)));
        A_vecs = readmatrix(allLines(i+(11:13)));
        A_vals = readvec(allLines{i+16}(9:end));
        A(n,:) = A_vals;
        if ~all(A_vecs==00,'all')
          AFrame(n,:) = eulang(A_vecs);
        else
          AFrame(n,:) = [0 0 0];
        end
        i = i + 18;
      end

    case '$ EPRNMR_QTensor'
      % colon is missing after "Number of stored nuclei"
      valuestr = regexp(allLines{idx0+4},'\d+$','match','once');
      nNuclei = str2double(valuestr);
      rho_present = contains(allLines{idx0+7},'true');
      d = 18;
      if rho_present
        d = 19;
      end
      i = idx0 + 8;
      Q = zeros(nNuclei,3);
      QFrame = zeros(nNuclei,3);
      for n = 1:nNuclei
        m = regexp(allLines{i},'(\d+)\W+(\w+)\W*$','tokens','once');
        atomno = str2double(m{1});
        element = m{2};
        Q_raw = readmatrix(allLines(i+(6:8)));
        Q_vecs = readmatrix(allLines(i+(11:13)));
        Q_vals = readvec(allLines{i+16}(12:end));
        Q(n,:) = Q_vals;
        QFrame(n,:) = eulang(Q_vecs);
        i = i + d;
      end

    otherwise
      fprintf('Skipping section %s\n',sectionTitle);

  end
end

% Read geometry
geometrySection = find(contains(allLines,'!GEOMETRY!'));
if ~isempty(geometrySection)
  idx0 = geometrySection;
  nAtoms = readvalue(allLines{idx0+1});
  idx = idx0+3;
  for iAtom = 1:nAtoms
    idx = idx + 1;
    xyz(iAtom,:) = readvec(allLines{idx}(22:end));
    atoms{iAtom} = regexp(allLines{idx}(1:22),'[a-zA-Z]+','match','once');
  end
end

% Remove nuclei with hyperfine coupling below cutoff
if ~isempty(A)
  keep = max(abs(A),[],2) >= abs(HyperfineCutoff);
  A = A(keep,:);
  AFrame = AFrame(keep,:);
  if ~isempty(Q)
    Q = Q(keep,:);
    QFrame = QFrame(keep,:);
  end
end

% Collect imported data into spin system structure
Sys = struct;
if ~isempty(S), Sys.S = S; end
if ~isempty(g), Sys.g = g; end
if ~isempty(gFrame), Sys.gFrame = gFrame; end
if ~isempty(D), Sys.D = D; end
if ~isempty(DFrame), Sys.DFrame = DFrame; end
if ~isempty(A), Sys.A = A; end
if ~isempty(AFrame), Sys.AFrame = AFrame; end
if ~isempty(Q), Sys.Q = Q; end
if ~isempty(QFrame), Sys.QFrame = QFrame; end
if ~isempty(E), Sys.Energy_Eh = E; end
if ~isempty(charge), Sys.charge = charge; end
if ~isempty(atoms), Sys.atoms = atoms; end
if ~isempty(xyz), Sys.xyz = xyz; end
end

%==========================================================================

function M = readmatrix(lines)
M(1,:) = sscanf(lines{1},'%g').';
M(2,:) = sscanf(lines{2},'%g').';
M(3,:) = sscanf(lines{3},'%g').';
M = M(:,2:end);
end

function val = readvalue(L)
i = find(L==':');
val = str2double(L(i+1:end));
end

function vec = readvec(lines)
if iscell(lines), lines = lines{1}; end
vec(1,:) = sscanf(lines,'%g').';
end
