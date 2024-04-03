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
    data(k).Element{iAtom} = tok{2};
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
  data(iStructure).Nucs = data(iStructure).Element;
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
      g_raw = readmatrix(allLines(idx0+(8:10)));
      g_sym = (g_raw.'*g_raw)^(1/2);
      recalcVecs = true;
      if recalcVecs
        g_sym = (g_sym+g_sym.')/2; % symmetrize numerically
        [g_vecs,g_vals] = eig(g_sym);
        g_vals = diag(g_vals).';
        if det(g_vecs)<0
          g_vecs(:,1) = -g_vecs(:,1);
        end
      else
        g_vecs = readmatrix(allLines(idx0+(13:15)));
        g_vals = readvec(allLines{idx0+18}(9:end));
      end
      data(iStructure).graw = g_raw;
      data(iStructure).g = g_vecs;
      data(iStructure).gvals = g_vals;
      data(iStructure).gFrame = eulang(g_vecs.');

    case '$ EPRNMR_DTensor'
      D_raw = readmatrix(allLines(idx0+(8:10)));
      recalcVecs = true;
      if recalcVecs
        [D_vecs,D_vals] = eig(D_raw);
        D_vals = diag(D_vals).';
        if det(D_vecs)<0
          D_vecs(:,1) = -D_vecs(:,1);
        end
      else
        D_vals = readvec(allLines{idx0+13}(9:end));
        D_vecs = readmatrix(allLines(idx0+(16:18)));
      end
      D_raw = D_raw*100*clight/1e6; % cm^-1 -> MHz
      D_vals = D_vals*100*clight/1e6; % cm^-1 -> MHz
      D_frame = eulang(D_vecs.');
      data(iStructure).Draw = D_raw;
      data(iStructure).Dvals = D_vals;
      data(iStructure).DFrame = D_frame;

    case '$ EPRNMR_ATensor'
      % colon is missing after "Number of stored nuclei"
      valuestr = regexp(allLines{idx0+4},'\d+$','match','once');
      nStoredNuclei = str2double(valuestr);
      A = cell(1,nAtoms);
      Araw = cell(1,nAtoms);
      AFrame = cell(1,nAtoms);
      i = idx0 + 7;
      for n = 1:nStoredNuclei
        m = regexp(allLines{i},'(\d+)\W+(\w+)\W*$','tokens','once');
        if isempty(m)
          break
        end
        atomno = str2double(m{1})+1; % 0-based -> 1 -based
        A_raw = readmatrix(allLines(i+(6:8)));
        A_vecs = readmatrix(allLines(i+(11:13)));
        A_vals = readvec(allLines{i+16}(9:end));
        Araw{atomno} = A_raw;
        A{atomno} = A_vals;
        if ~all(A_vecs(:)==0)
          AFrame{atomno} = eulang(A_vecs.');
        else
          AFrame{atomno} = [0 0 0];
        end
        i = i + 18;
      end
      data(iStructure).Araw = Araw;
      data(iStructure).A = A;
      data(iStructure).AFrame = AFrame;

    case '$ EPRNMR_QTensor'
      warning('EPRNMR_QTensor is not supported for property.txt files. Use the main output file instead.');

    case '$ EPRNMR_EFGTensor'
      valuestr = regexp(allLines{idx0+4},'\d+$','match','once');
      nStoredNuclei = str2double(valuestr);
      rho_present = contains(allLines{idx0+7},'true');
      d = 18;
      if rho_present
        d = 19;
      end
      i = idx0 + 8;
      Q = cell(1,nAtoms);
      Qraw = cell(1,nAtoms);
      QFrame = cell(1,nAtoms);
      for n = 1:nStoredNuclei
        m = regexp(allLines{i},'(\d+)\W+(\w+)\W*$','tokens','once');
        atomno = str2double(m{1})+1;  % 0-based -> 1-based
        Q_raw = readmatrix(allLines(i+(6:8)));
        Q_vecs = readmatrix(allLines(i+(11:13)));
        Q_vals = readvec(allLines{i+16}(12:end));
        Qraw{atomno} = Q_raw;
        Q{atomno} = Q_vals;
        QFrame{atomno} = eulang(Q_vecs.');
        i = i + d;
      end
      data(iStructure).Q = Q;
      data(iStructure).Qraw = Qraw;
      data(iStructure).QFrame = QFrame;

    otherwise
      %fprintf('Skipping section %s\n',sectionTitle);

  end
end

if ~isfield(data,'S')
  error('Spin multiplicity is missing in property file.');
end

% Copy relevant data to spin system structure
%--------------------------------------------------------------------------
for iStructure = nStructures:-1:1
  d = data(iStructure);
  
  % Coordinates
  if ~isempty(d.xyz)
    Sys(iStructure).xyz = d.xyz;
  end

  % Spin multiplicity
  Sys(iStructure).S = d.S;

  % g tensor
  if isfield(d,'g')
    if ~isempty(d.g)
      Sys(iStructure).g = d.gvals;
    end
    if ~isempty(d.gFrame)
      Sys(iStructure).gFrame = d.gFrame;
    end
  end

  % D tensor
  if isfield(d,'Dvals')
    if ~isempty(d.Dvals)
      Sys(iStructure).D = d.Dvals;
    end
    if ~isempty(d.DFrame)
      Sys(iStructure).DFrame = d.DFrame;
    end
  end

  % Compile nuclear data (isotopes, hyperfine coupling, quuadropole coupling)
  idx = 0;
  for iAtom = 1:nAtoms
    Agiven = isfield(d,'A') && ~isempty(d.A{iAtom});
    Qgiven = isfield(d,'Q') && ~isempty(d.Q{iAtom});
    if Agiven || Qgiven
      idx = idx + 1;
      Sys(iStructure).Nucs{idx} = d.Element{iAtom};
      Sys(iStructure).NucsIdx(idx) = iAtom;
      if Agiven
        Sys(iStructure).A(idx,1:3) = d.A{iAtom};
        Sys(iStructure).AFrame(idx,1:3) = d.AFrame{iAtom};
      end
      if Qgiven
        Sys(iStructure).Q(idx,1:3) = d.Q{iAtom};
        Sys(iStructure).QFrame(idx,1:3) = d.QFrame{iAtom};
      end
    end
  end

  % Store all other data in spin system structure
  Sys(iStructure).data = data;
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
