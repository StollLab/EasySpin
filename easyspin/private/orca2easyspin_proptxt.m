% orca2easyspin_maintxt   Read EPR properties from main ORCA output file

function Sys = orca2easyspin_proptxt(propfilename)

% Read entire file into cell array
fh = fopen(propfilename);
allLines = textscan(fh,'%s','whitespace','','delimiter','\n');
allLines = allLines{1};
fclose(fh);

if ~contains(allLines{2},"!PROPERTIES!")
  error('This is not a valid text-based ORCA property file.');
end

% Read geometry
%--------------------------------------------------------------------------
geometrySection = find(contains(allLines,'!GEOMETRIES!'),1);
if isempty(geometrySection)
  error('Cannot find !GEOMETRIES! section.')
end
geometrySection = find(contains(allLines,'!GEOMETRY!'));
nStructures = numel(geometrySection);
for k = 1:numel(geometrySection)
  idx0 = geometrySection(k);
  nAtoms = readvalue(allLines{idx0+1});
  data(k).nAtoms = nAtoms;
  idx = idx0+3;
  for iAtom = 1:nAtoms
    idx = idx + 1;
    data(k).xyz(iAtom,:) = readvec(allLines{idx}(22:end));
    tok = regexp(allLines{idx}(1:22),'(\d+)\s+([a-zA-Z]+)','tokens','once');
    data(k).atoms{iAtom} = tok{2};
  end
end


% Loop over all sections and parse relevant information
%--------------------------------------------------------------------------
% Locate lines with section titles
sections = find(cellfun(@(L)L(1)=='$',allLines));
for iSection = 1:numel(sections)
  idx0 = sections(iSection);
  sectionTitle = allLines{idx0};
  iStructure = str2double(regexp(allLines{sections(iSection)+2},'(\d+)$','match','once'));
  data(iStructure).Nucs = data(iStructure).atoms;
  switch sectionTitle
    case '$ SCF_Energy'
      data(iStructure).E = readvalue(allLines{idx0+4});

    case '$ DFT_Energy'
      % not relevant

    case '$ Mayer_Pop'
      % not relevant

    case '$ Calculation_Info'
      multiplicity = readvalue(allLines{idx0+4});
      data(iStructure).S = (multiplicity-1)/2;
      data(iStructure).charge = readvalue(allLines{idx0+5});

    case '$ SCF_Electric_Properties'
      % not relevant

    case '$ EPRNMR_GTensor'
      %g_raw = readmatrix(allLines(idx0+(8:10)));
      g_vecs = readmatrix(allLines(idx0+(13:15)));
      g_vals = readvec(allLines{idx0+18}(9:end));
      data(iStructure).g = g_vals;
      data(iStructure).gFrame = eulang(g_vecs);

    case '$ EPRNMR_DTensor'
      D_raw = readmatrix(allLines(idx0+(8:10)));
      D_vals = readvec(allLines{idx0+13}(9:end));
      D_vecs = readmatrix(allLines(idx0+(16:18)));
      %[D_vecs,D_vals] = eig(D_raw);
      %D_vals = diag(D_vals).';
      data(iStructure).D = D_vals;
      data(iStructure).DFrame = eulang(D_vecs);

    case '$ EPRNMR_ATensor'
      % colon is missing after "Number of stored nuclei"
      valuestr = regexp(allLines{idx0+4},'\d+$','match','once');
      nNuclei = str2double(valuestr);
      A = zeros(nNuclei,3);
      AFrame = zeros(nNuclei,3);
      i = idx0 + 7;
      for n = 1:nNuclei
        m = regexp(allLines{i},'(\d+)\W+(\w+)\W*$','tokens','once');
        if isempty(m)
          break
        end
        atomno = str2double(m{1}); % 0-based
        element = m{2};
        Nucs{atomno+1} = element;
        %A_raw = readmatrix(allLines(i+(6:8)));
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
      data(iStructure).A = A;
      data(iStructure).AFrame = AFrame;

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
        Nucs{atomno+1} = element;
        %Q_raw = readmatrix(allLines(i+(6:8)));
        Q_vecs = readmatrix(allLines(i+(11:13)));
        Q_vals = readvec(allLines{i+16}(12:end));
        Q(n,:) = Q_vals;
        QFrame(n,:) = eulang(Q_vecs);
        i = i + d;
      end
      data(iStructure).Q = Q;
      data(iStructure).QFrame = QFrame;

    otherwise
      %fprintf('Skipping section %s\n',sectionTitle);

  end
end

% Collect imported data into spin system structure

Sys = data;
for iStructure = 1:nStructures
  Sys(iStructure).Nucs = data(iStructure).atoms;
end

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
