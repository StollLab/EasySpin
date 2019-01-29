% curry  Computation of magnetometry data
%         (magnetic moment, susceptibility)
%
%   curry(Sys,Exp)
%   curry(Sys,Exp,Opt)
%   muz = curry(...)
%   [muz,chizz] = curry(...)
%
%    Calculates the magnetic moment and the molar static magnetic
%    susceptibility of a powder sample for given values of
%    applied magnetic field and temperature.
%
%    Input:
%      Sys    spin system
%      Exp    experimental parameter
%        Field        list of field values (mT)
%        Temperature  list of temperatures (K)
%        
%        mField, mTemperature 
%        chiField, chiTemperature
%                     list of field values and temperatures for the
%                     calculation of the magnetic moment (mField &
%                     mTemperature) and the magnetic susceptibility
%                     (chiField & chiTemperature)
%
%
%      Opt    calculation options
%        nKnots       number of knots for powder average
%        Output       list of keywords defining the output and the order of
%                     it
%                     the following keywords are allowed:
%                     'MvsB', 'MvsBCGS', 'MvsBSI', 'Chi',
%                     'ChiT', '1overChi', 'MuEff','ChiSI'
%                     'ChiTSI', '1overChiSI', 'MuEffSI',
%                     'ChiCGS', 'ChiTCGS', '1overChiCGS', 
%                     'MuEffCGS', 'OneColoumn'
%        Method       'operator' (default) or 'energies'
%                     calculation method
%  
%                     
%    Output (if Opt.Output is not given):
%      muz      magnetic moment along zL axis
%               in units of Bohr magnetons
%      chizz    molar susceptibility, zLzL component
%               in SI units (m^3 mol^-1)
%    
%    zL is the direction of the applied static magnetic field
%
%    The size of muz and chizz is nB x nT, where nB is the number of
%    field values in Exp.Field and nT is the number of temperature values
%    in Exp.Temperature.
%
%   If no output argument is given, the computed data are plotted.

function varargout = curry(Sys,Exp,Opt)

if (nargin==0), help(mfilename); return; end

if (nargin<2), Exp = struct; end
if (nargin<3), Opt = struct; end
if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;

logmsg(1,['=begin=curry======' datestr(now) '=================']);

calculateChi = (nargout==0) || (nargout>1);
doPlot = (nargout==0);

% Spin system
%-------------------------------------------------
logmsg(1,'Spin system');
if iscell(Sys)
  error('curry does not support calculations with multiple components.');
end
if ~isfield(Sys,'TIP'), Sys.TIP = 0; end

% Experimental parameters
%-------------------------------------------------
logmsg(1,'Experimental parameters');

% different parameters for chi and m
cmpstr = {'chiTemperature', 'chiField', 'mTemperature', 'mField'};
check = isfield(Exp,cmpstr);
diffParam = any(check);
if ~diffParam
  % Field
  if ~isfield(Exp,'Field')
    Exp.Field = 0;
    disp('Exp.Field is missing, assuming zero field.');
  end
  B = Exp.Field/1e3; % magnetic field, T
  nFields = numel(B);
  logmsg(1,'  number of field values: %d',nFields);
  
  % Temperature
  if ~isfield(Exp,'Temperature')
    error('Exp.Temperature is missing.');
  end
  
  T = Exp.Temperature(:).'; % temperature, K
  nTemperatures = numel(T);
  logmsg(1,'  number of temperature values: %d',nTemperatures);
  
  if any(T<0)
    error('Negative temperatures are not possible.')
  end
  zeroTemp = T==0;
end

% Other fields
if ~isfield(Exp,'CrystalOrientation')
  Exp.CrystalOrientation = [];
end
if ~isfield(Exp,'CrystalSymmetry')
  Exp.CrystalSymmetry = '';
end
if ~isfield(Exp,'MolFrame')
  Exp.MolFrame = [];
end

PowderSimulation = isempty(Exp.CrystalOrientation);

if PowderSimulation
  logmsg(1,'Powder calculation');
else
  logmsg(1,'Crystal calculation');
end

% Options
%-------------------------------------------------
logmsg(1,'Options');
if ~isfield(Opt,'nKnots')
  Opt.nKnots = 10;
end
logmsg(1,'  number of knots: %d',Opt.nKnots);
if ~isfield(Opt,'Symmetry')
  Opt.Symmetry = []; % needed for p_symandgrid
end
if ~isfield(Opt,'SymmFrame')
  Opt.SymmFrame = []; % needed for p_symandgrid
end

calculateM = true;
OneColoumn = false;
if ~isfield(Opt,'Output')
  if nargout == 2
    out = logical([1,zeros(1,14);zeros(1,3),1,zeros(1,11)]);
    len =2;
  elseif nargout == 1
    out = logical([1,zeros(1,14)]);
    len = 1;
  end
else
  keywords = {'MvsB', 'MvsBCGS', 'MvsBSI' ...
    'Chi', 'ChiT', '1overChi', 'MuEff', ...
    'ChiSI', 'ChiTSI', '1overChiSI', 'MuEffSI', ...
    'ChiCGS', 'ChiTCGS', '1overChiCGS', 'MuEffCGS', ...
    'OneColoumn'};
  r = textscan(Opt.Output,'%s ');
  len = length(r{1});
  for n = len:-1:1
      out(n,:) =strcmpi(r{1}{n},keywords);  
  end
  if~(all(any(out,2)))
    error([r{1}{~any(out,2)}, ' is not an allowed keyword for Opt.Output!'])
  end
  for k = find(out(:,16))
    if isempty(k), continue; end;
    OneColoumn = true;
    out = out(~(1:len==k),:); % remove OneColumn from Outputlist
    len = len -1;
    logmsg(1,'  single coloumn output');
    if isempty(out), error('Specify which output should appear as one Coloumn!'); end
  end
  out = out(:,1:15); % remove OneColumn from Outputlist
  if ~OneColoumn && ~(nargout==len)
    error('Number of output variables do not match requested outputs!');
  end
  logmsg(1,'  number of outputs: %d',len);  
  if any(any(out(:,4:15))), calculateChi = true; end
  if ~any(any(out(:,1:3))), calculateM = false; end 
end

if isfield(Opt,'Method')
  if strcmpi(Opt.Method,'operator')
    DoOp = true;
  elseif strcmpi(Opt.Method,'energies')
    DoOp = false;
    if calculateChi && ~diffParam && any(B(:))
      error('Energie method can calculate susceptibility only in the absence of magnetic field!');
    end
  else
    error('Opt.Method can be Operator or Energies!');
  end
else
  DoOp = true;
  Opt.Method = 'Operator';
end
logmsg(1,'  calculation method: %s',Opt.Method);

%----- differnt fields and temperatures for m and chi
if diffParam
  doCalc = false;
  if isfield(Opt,'Output')
    OptS = rmfield(Opt,'Output');
  else
    OptS = Opt;
  end
  if calculateM
    if all(check(3:4))
      ExpM = Exp;
      ExpM.Temperature = Exp.mTemperature;
      ExpM.Field = Exp.mField;
      ExpM = rmfield(ExpM,cmpstr);
      muz = curry(Sys,ExpM,OptS)*bmagn;
    else
      error('mTemperature and mField have to be given!');
    end
  end
  if calculateChi
    if all(check(1:2))
      ExpChi = Exp;
      ExpChi.Temperature = Exp.chiTemperature;
      ExpChi.Field = Exp.chiField;
      ExpChi = rmfield(ExpChi,cmpstr);
      [t, chizz_SI] = curry(Sys,ExpChi,OptS);
      clear t;
    else
      error('chiTemperature and chiField have to be given!');
    end
  end
else
    doCalc = true;
end
if(doCalc)
  % Set up Hamiltonian and magnetic dipole moment
  %-------------------------------------------------
  % zero-field Hamiltonian F (MHz)
  % magnetic dipole moment operators muOpxM,muOpyM,muOpzM (MHz/mT)
  % all are in the molecular frame
  [H0,GxM,GyM,GzM] = sham(Sys);
  
  % zero-field spin Hamiltonian
  H0 = H0*1e6*planck; % MHz -> J
  
  % magnetic dipole operators, in molecular frame
  muOpxM = -GxM*1e6*1e3*planck; % MHz/mT -> J/T
  muOpyM = -GyM*1e6*1e3*planck; % MHz/mT -> J/T
  muOpzM = -GzM*1e6*1e3*planck; % MHz/mT -> J/T
  
  % Set up sample orientations
  %-------------------------------------------------
  Exp.PowderSimulation = PowderSimulation; % for communication with p_*
  [Exp,Opt] = p_symandgrid(Sys,Exp,Opt);
  % Process crystal orientations, crystal symmetry, and frame transforms
  % This sets Orientations, nOrientations, nSites and AverageOverChi
  [Orientations,nOrientations,nSites,AverageOverChi] = p_crystalorientations(Exp,Opt);
  Exp.OriWeights = Exp.OriWeights/4/pi;
  
  beta = 1./T/boltzm;
  
  % Initialize output arrays
  muz = zeros(nFields,nTemperatures);
  if calculateChi
    chizz = zeros(nFields,nTemperatures);
  end
  
  % Calculation loop
  %-------------------------------------------------
  % all calculations are done in the molecular frame
  if DoOp
    for iOri = 1:nOrientations
      [xLab_M,yLab_M,zLab_M] = erot(Orientations(iOri,:),'rows');
      % projection of magnetic moment operator onto lab axes
      % (field direction along z axis of lab frame, zLab)
      muOpzL = zLab_M(1)*muOpxM + zLab_M(2)*muOpyM + zLab_M(3)*muOpzM; % J/T
      
      for iB = 1:nFields
        
        [V,E] = eig(H0 - B(iB)*muOpzL);
        E = diag(E) - E(1); % J
        populations = exp(-E*beta);
        if zeroTemp, populations(1) = 1; end
        
        % expectation values for each spin state
        muz_expect = real(diag(V'*muOpzL*V));
        
        % population-weighted average
        muz_avg = (muz_expect.'*populations)./sum(populations,1);
        
        % accumulation for powder average
        muz(iB,:) = muz(iB,:) + Exp.OriWeights(iOri)*muz_avg;
        
        if calculateChi
          dB = eps^(1/3)*max(B(iB),1); % optimal step size for numerical derivative
          B1 = B(iB) + dB; dB = B1 - B(iB); % prevent round-off errors
          
          [V,E] = eig(H0 - (B(iB)+dB)*muOpzL);
          E = diag(E) - E(1); % J
          populations = exp(-E*beta);
          if zeroTemp, populations(1) = 1; end
          
          % population-weighted average
          muz_expect = real(diag(V'*muOpzL*V));
          muz_avg2 = (muz_expect.'*populations)./sum(populations,1);
          
          chizz_ = (muz_avg2-muz_avg)/dB;
          chizz(iB,:) = chizz(iB,:) + Exp.OriWeights(iOri)*chizz_;
        end
        
      end
    end
  else
    if calculateChi
      E0 = eig(H0);
      mEchi = E0(1);
      E0 = E0- mEchi;
      lE = length(E0);
      betachi = repmat(beta,lE,1);
      lZ0 = log(sum(exp(-repmat(E0,1,nTemperatures).*betachi)));
    end
    for iOri = 1:nOrientations
      %calculate M and chi with approximations
      [xLab_M,yLab_M,zLab_M] = erot(Orientations(iOri,:),'rows');
      muOpzL = zLab_M(1)*muOpxM + zLab_M(2)*muOpyM + zLab_M(3)*muOpzM;
      % projection of magnetic moment operator onto lab axes
      % (field direction along z axis of lab frame, zLab)
      if calculateM
        for iB = nFields:-1:1
          dB = eps^(1/3)*max(B(iB),1); % optimal step size for numerical derivative
          B1 = B(iB) + dB; dB = B1 - B(iB); % prevent round-off errors
          B2 = B(iB) - dB;
          E1 = eig(H0+ B1*muOpzL);
          E2 = eig(H0+ B2*muOpzL);
          mE = min([E1;E2]);
          E1 = E1-mE;
          E2 = E2-mE;
          lZ1 = log(sum(exp(-E1*beta)));
          lZ2 = log(sum(exp(-E2*beta)));
          mz(:,iB)= boltzm * Exp.Temperature .* (lZ1-lZ2).'/(2*dB);
        end
        % accumulation for powder average
        muz = muz + Exp.OriWeights(iOri)*mz.';
      end
      if calculateChi
        dB = eps^(1/3)*1e3;
        E1 = eig(H0+ dB*muOpzL);
        E1 = E1-mEchi;
        lZ1 = log(sum(exp(-repmat(E1,1,nTemperatures).*betachi)));
        chizz_ = boltzm *Exp.Temperature.' .* (lZ1-lZ0)*2/dB^2;
        chizz = chizz + Exp.OriWeights(iOri)*chizz_;
      end
    end
  end
  
  % Unit conversions
  if  calculateChi
    chizz_SI = chizz*mu0*avogadro +Sys.TIP;   % single molecule SI -> molar SI
  end
  
  %chizz_cgs = chizz_SI/(4*pi*1e-6);   % SI -> CGS-emu unit conversion
  
  % Graphical plotting
  %-----------------------------------------------------
  if doPlot
    DataDimensionality = (nFields~=1) + (nTemperatures~=1);
    if DataDimensionality==0
      % do nothing
      clf
      
    elseif DataDimensionality==1
      if (nFields==1)
        x = T;
        xLab = 'temperature (K)';
      else
        x = B;
        xLab = 'magnetic field (T)';
      end
      
      cla
      subplot(2,2,1)
      plot(x,muz/bmagn);
      axis tight
      ylabel('\mu_z (\mu_B), \mu_{mol,z} (N_A\mu_B)')
      
      subplot(2,2,2)
      plot(x,muz*avogadro);
      axis tight
      ylabel('\mu_{mol,z} (J T^{-1} mol^{-1})')
      
      subplot(2,2,3)
      plot(x,chizz_SI);
      ylabel('\chi_{mol,zz} (m^3 mol^{-1})');
      
      subplot(2,2,4)
      plot(x,chizz_SI.*T);
      ylabel('\chi_{mol,zz}T (K m^3 mol^{-1})');
      
      for i=1:4
        subplot(2,2,i);
        xlabel(xLab);
        xlim([min(x) max(x)]);
      end
      
    else
      cla
      Plot2D = (nFields>10);
      chizzT_SI = chizz_SI.*repmat(T(:).',nFields,1);
      if Plot2D
        subplot(2,2,1);
        surf(T,B,muz/bmagn);
        shading flat
        zlabel(' \mu_z (\mu_B), \mu_{mol,z} (N_A\mu_B)')
        subplot(2,2,2);
        surf(T,B,muz*avogadro);
        shading flat
        zlabel(' \mu_{mol,z} (J T^{-1} mol^{-1})')
        subplot(2,2,3);
        surf(T,B,chizz_SI);
        shading flat
        zlabel('\chi_{mol,zz} (m^3 mol^{-1})');
        subplot(2,2,4);
        surf(T,B,chizzT_SI);
        shading flat
        zlabel('\chi_{mol,zz}T (K m^3 mol^{-1})');
        for i=1:4
          subplot(2,2,i);
          xlabel('T (K)');
          ylabel('B (T)')
          xlim([min(T) max(T)]);
          ylim([min(B) max(B)]);
        end
      else
        subplot(2,2,1)
        plot(T,muz/bmagn);
        axis tight
        ylabel('\mu_z (\mu_B), \mu_{mol,z} (N_A\mu_B)')
        
        subplot(2,2,2)
        plot(T,muz*avogadro);
        axis tight
        ylabel('\mu_{mol,z} (J T^{-1}mol^{-1})')
        
        subplot(2,2,3)
        plot(T,chizz_SI);
        ylabel('\chi_{mol,zz} (m^3 mol^{-1})');
        
        subplot(2,2,4)
        plot(T,chizzT_SI);
        ylabel('\chi_{mol,zz}T (K m^3 mol^{-1})');
        
        for i=1:4
          subplot(2,2,i);
          xlabel('temperature (K)');
          xlim([min(T) max(T)]);
        end
        
      end
    end
  end
end
% Assign output arguments
%---------------------------------------------------------------
if OneColoumn
  if diffParam
    nFields = numel(Exp.chiField);
    T = Exp.chiTemperature;
  else
    outdim =  nFields*nTemperatures;
  end
  for n =len:-1:1
    switch find(out(n,:))
      case {1,2}, temp = muz/bmagn; %MvsB,MvsBCGS
      case 3, temp = muz *avogadro; %MvsBSi
      case {4,8}, temp = chizz_SI; %Chi,ChiSI
      case {5,9}, temp = chizz_SI.*repmat(T(:).',nFields,1); %ChiT,ChiTSi
      case {6,10}, temp = 1./chizz_SI; %1overChi,1overChiSI
      case {7,11}, temp = sqrt(chizz_SI.*repmat(T(:).',nFields,1)/8); %MuEff, MuEffSI
      case 12, temp = chizz_SI/(4*pi*1e-6); %ChiCGS
      case 13, temp = chizz_SI.*repmat(T(:).',nFields,1)/(4*pi*1e-6); %ChiTCGS
      case 14, temp = 1./chizz_SI*(4*pi*1e-6); %1overChiCGS
      case 15, temp = sqrt(chizz_SI.*repmat(T(:).',nFields,1)/8/(4*pi*1e-6)); %MuEffCGS
    end
    if diffParam
      tcell{n} = temp(:);
    else
      for m= nTemperatures:-1:1
        templine((m-1)*nFields+1:m*nFields) = temp(:,m);
      end
      varargout{1}((n-1)*outdim+1:n*outdim) = templine;
    end
  end
  if diffParam
     varargout{1} =[];
    for m=1:length(tcell)
      varargout{1} = [varargout{1}, tcell{m}.'];
    end
  end
else
  for n =len:-1:1
    % different exp. parameter for m and chi, adjust nFields accordingly
    if diffParam
      nFields = numel(Exp.chiField); 
      T = Exp.chiTemperature;
    end
    switch find(out(n,:))
      case {1,2}, varargout{n} = muz/bmagn; %MvsB,MvsBCGS
      case 3, varargout{n} = muz *avogadro; %MvsBSi
      case {4,8}, varargout{n} = chizz_SI; %Chi,ChiSI
      case {5,9}, varargout{n} = chizz_SI.*repmat(T(:).',nFields,1); %ChiT,ChiTSi
      case {6,10}, varargout{n} = 1./chizz_SI; %1overChi,1overChiSI
      case {7,11}, varargout{n} = sqrt(chizz_SI.*repmat(T(:).',nFields,1)*8); %MuEff, MuEffSI
      case 12, varargout{n} = chizz_SI/(4*pi*1e-6); %ChiCGS
      case 13, varargout{n} = chizz_SI.*repmat(T(:).',nFields,1)/(4*pi*1e-6); %ChiTCGS
      case 14, varargout{n} = 1./chizz_SI*(4*pi*1e-6); %1overChiCGS
      case 15, varargout{n} = sqrt(chizz_SI.*repmat(T(:).',nFields,1)*8/(4*pi*1e-6)); %MuEffCGS
    end
  end
end

% switch (nargout)
%   case 0
%   case 1
%     varargout = {muz/bmagn};
%   case 2
%     varargout = {muz/bmagn,chizz_SI};
%   case 3
%     varargout = {muz/bmagn,chizz_SI,chizzT_SI};
%   otherwise
%     varargout = cell(1,nargout);
%     varargout(1:3) = {muz/bmagn,chizz_SI,chizzT_SI};
% end

logmsg(1,'=end=curry========%s=================\n',datestr(now));

clear global EasySpinLogLevel
