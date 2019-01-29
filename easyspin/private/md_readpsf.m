%  md_readpsf  Read a PSF (protein structure file, in CHARMM or 
%                    X-PLOR format) to obtain nitroxide spin label 
%                    information.
%
%   psf = md_readpsf(FileName, ResName, AtomNames);
%
%   Input:
%     FileName       character array
%                    Name of DCD file to be read.
%
%     SegName        character array
%                    Name of segment consisting of the protein and spin
%                    label.
%
%     ResName        character array
%                    Name of residue assigned to spin label side chain,
%                    e.g. "CYR1" is the default used by CHARMM-GUI.
%
%     AtomNames      structure array
%                    Structure array containing the atom names used in the 
%                    PSF to refer to the following atoms in the nitroxide 
%                    molecule:
%
%                                    ON (ONname)
%                                    |
%                                    NN (NNname)
%                                  /   \
%                        (C1name) C1    C2 (C2name)
%                                 |     |
%                       (C1Rname) C1R = C2R (C2Rname)
%                                 |
%                       (C1Lname) C1L
%                                 |
%                       (S1Lname) S1L
%                                /
%                      (SGname) SG
%                               |
%                      (CBname) CB
%                               |
%                   (Nname) N - CA (CAname)
%
%   Output:
%     psf            structure array
%                    NATOM: number of atoms to be included in simulation
%                    idx_Protein: indices of the protein and label atoms
%                    idx_ProteinCA: indices of protein's alpha carbon atoms
%                    idx_SpinLabel: indices of the spin label's atoms
%                    idx_ON: index corresponding to atom ONname
%                    idx_NN: " " NNname
%                    idx_C1: " " C1name
%                    idx_C2: " " C2name
%                    idx_C1R: " " C1Rname
%                    idx_C2R: " " C2Rname
%                    idx_C1L: " " C1Lname
%                    idx_S1L: " " S1Lname
%                    idx_SG: " " SGname
%                    idx_CB: " " CBname
%                    idx_CA: " " CAname
%                    idx_N: " " Nname
%
% For more information on PSF file structure, see:
%  http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html

function psf = md_readpsf(FileName, SegName, ResName, AtomNames)

% Note:
% This function assumes the following order of sections in the PSF, which 
% is standard:
%   !NTITLE
%   !NATOM
%   !NBOND
%   ...

if ~ischar(FileName)
  error('FileName must be given as a character array.')
end

if ~strcmpi(FileName(end-3:end),'.psf')
  error('Please give the full filename including the ".psf" extension.')
end

if ~ischar(ResName)
  error('resName must be given as a character array.')
end


% check AtomNames

if any(~cellfun(@ischar,fieldnames(AtomNames)))
  error('All fields in AtomNames must be given as character arrays.')
end

if ~isfield(AtomNames,'ONname')||~isfield(AtomNames,'NNname')...
    ||~isfield(AtomNames,'C1name')||~isfield(AtomNames,'C2name')
  error('One or more required atom names is missing. Please check documentation.')
end

% open the file
FileID = fopen(FileName, 'r');
if FileID<1, error('File "%s" could not be opened.', FileName); end
ensurefclose = onCleanup(@() fclose(FileID));

% proper PSF formatting dictates that the first line of the file should
% start with "PSF"
line = fgetl(FileID);
if ~ischar(line)||~strcmpi(line(1:3),'PSF')
  error('File "%s" does not have proper PSF format. See documentation for details.', FileName)
end

% initialization
% -------------------------------------------------------------------------

% psf.isCHARMM = 0;  % CHARMM and X-PLOR formats only differ in how dihedral
%                    % force field terms are listed, which is irrelevant for
%                    % reading atom data

psf.NATOM = 0;  % this parameter should be checked against the trajectory 
                % file to ensure that the files are compatible
                
AtomFormat = '%10d %8s %8d %8s %8s %4s %14f%14f%8d %*f %*f';

reachedNTITLE = 0;
reachedNATOM = 0;

% read the file
% -------------------------------------------------------------------------

while ~feof(FileID)
  line = strtrim(fgetl(FileID));
  
  if regexp(line, '.*\d.*!\w.*')
    
    line = strsplit(line);
    
    % each section should begin with a line in the following format:
    %    [nLines] ![section]...
    % where nLines gives the number of lines that the section occupies 
    % after this beginning line
    nLines = round(str2double(line{1}));
    section = line{2};
    
    if ~isempty(strfind(section,'NTITLE'))
      % skip this section
      reachedNTITLE = 1;
    end
    
    if ~isempty(strfind(section,'NATOM'))
      if ~reachedNTITLE
        error('Section ordering in "%s" is not standard. See documentation for proper formatting.', FileName)
      end
      reachedNATOM = 1;
      psf.NATOM = nLines;
      
      % read NATOM section
      FileContents = textscan(FileID, AtomFormat, nLines);
      segmentNames = FileContents(2);
      segmentNames = segmentNames{1};
      residueNames = FileContents(4);
      residueNames = residueNames{1};
      atomNames = FileContents(5);
      atomNames = atomNames{1};
      mass = FileContents(8);
      psf.mass = mass{1};
      
      % filter for atoms belonging to the protein and spin label
      idx_ProteinLabel = strcmpi(segmentNames,SegName);
      segmentNames = segmentNames(idx_ProteinLabel);
      residueNames = residueNames(idx_ProteinLabel);
      atomNames = atomNames(idx_ProteinLabel);
      
      % generate logical arrays to filter for protein alpha carbon and spin
      % label atoms
      psf.idx_ProteinCA = strcmpi(segmentNames,SegName) & strcmpi(atomNames,'CA');
      
      idx_SpinLabel = strcmpi(segmentNames,SegName) & strncmpi(residueNames,ResName,4);
      psf.idx_ON = strcmpi(atomNames(idx_SpinLabel),AtomNames.ONname);
      psf.idx_NN = strcmpi(atomNames(idx_SpinLabel),AtomNames.NNname);
      psf.idx_C1 = strcmpi(atomNames(idx_SpinLabel),AtomNames.C1name);
      psf.idx_C2 = strcmpi(atomNames(idx_SpinLabel),AtomNames.C2name);
      psf.idx_C1R = strcmpi(atomNames(idx_SpinLabel),AtomNames.C1Rname);
      psf.idx_C2R = strcmpi(atomNames(idx_SpinLabel),AtomNames.C2Rname);
      psf.idx_C1L = strcmpi(atomNames(idx_SpinLabel),AtomNames.C1Lname);
      psf.idx_S1L = strcmpi(atomNames(idx_SpinLabel),AtomNames.S1Lname);
      psf.idx_SG = strcmpi(atomNames(idx_SpinLabel),AtomNames.SGname);
      psf.idx_CB = strcmpi(atomNames(idx_SpinLabel),AtomNames.CBname);
      psf.idx_CA = strcmpi(atomNames(idx_SpinLabel),AtomNames.CAname);
      psf.idx_N = strcmpi(atomNames(idx_SpinLabel),AtomNames.Nname);
      psf.idx_SpinLabel = idx_SpinLabel;
      
      % the ProteinLabel indices will be read directly from the binary 
      % files, so convert to integer indices
      psf.idx_ProteinLabel = find(idx_ProteinLabel);
      
%       psf.idx_ProteinCA = find(idx_ProteinCA);
% %       psf.idx_ProteinCA = find(idx_ProteinCA) - psf.idx_ProteinLabel(1) + 1;
%       psf.idx_SpinLabel = find(idx_SpinLabel);
% %       psf.idx_SpinLabel = find(idx_SpinLabel) - psf.idx_ProteinLabel(1) + 1;
% 
%       psf.idx_ON = find(idx_ON);
%       psf.idx_NN = find(idx_NN);
%       psf.idx_C1 = find(idx_C1);
%       psf.idx_C2 = find(idx_C2);
%       psf.idx_C1R = find(idx_C1R);
%       psf.idx_C2R = find(idx_C2R);
%       psf.idx_C1L = find(idx_C1L);
%       psf.idx_S1L = find(idx_S1L);
%       psf.idx_SG = find(idx_SG);
%       psf.idx_CB = find(idx_CB);
%       psf.idx_CA = find(idx_CA);
%       psf.idx_N = find(idx_N);
      
%       psf.idx_ON = nonzeros(psf.idx_SpinLabel .* idx_ON);
%       psf.idx_NN = nonzeros(psf.idx_SpinLabel .* idx_NN);
%       psf.idx_C1 = nonzeros(psf.idx_SpinLabel .* idx_C1);
%       psf.idx_C2 = nonzeros(psf.idx_SpinLabel .* idx_C2);
%       psf.idx_C1R = nonzeros(psf.idx_SpinLabel .* idx_C1R);
%       psf.idx_C2R = nonzeros(psf.idx_SpinLabel .* idx_C2R);
%       psf.idx_C1L = nonzeros(psf.idx_SpinLabel .* idx_C1L);
%       psf.idx_S1L = nonzeros(psf.idx_SpinLabel .* idx_S1L);
%       psf.idx_SG = nonzeros(psf.idx_SpinLabel .* idx_SG);
%       psf.idx_CB = nonzeros(psf.idx_SpinLabel .* idx_CB);
%       psf.idx_CA = nonzeros(psf.idx_SpinLabel .* idx_CA);
%       psf.idx_N = nonzeros(psf.idx_SpinLabel .* idx_N);
    end
    
    if ~isempty(strfind(section,'NBOND'))
      if ~reachedNTITLE||~reachedNATOM
        error('Section ordering in "%s" is not standard. See documentation for proper formatting', FileName)
      end
      % we are done
      break
    end
  
  end
end

end