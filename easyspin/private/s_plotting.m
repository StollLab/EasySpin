function s_plotting(TimeAxis,Signal,Exp,Opt)
% This function plots the output of saffron and spidyan

if ~isfield(Exp,'DetOperator')
  % this field will be used as title in the plotting
  Exp.DetOperator = {'electron coherence'};
end

nDetOps = length(Exp.DetOperator);

if isfield(Exp,'nPoints')
  nDataPoints = prod(Exp.nPoints);
else
  nDataPoints = 1;
  if  ~Opt.SinglePointDetection && size(Signal,1) == length(Exp.DetOperator)
    % make first dimension the transient - required for the plotting
    % further down
    Signal = Signal.';
  end
end

% Set up figures
logmsg(1,'  setting up figures...');
% only make figures if a) transients were detected or b) single point
% detection with indirect dimension (but not more than 2) and more
% than one acquisition point
if ~Opt.SinglePointDetection || (Opt.SinglePointDetection && isfield(Exp,'nPoints') && length(Exp.nPoints) < 3 && nDataPoints > 1)
  LabelsDetectionOp = cell(1,nDetOps);
  for iDetOp = 1 : nDetOps
    figure
    if ischar(Exp.DetOperator{iDetOp})
      % title for detection operators that used the sop syntax
      title(['Signal of ' Exp.DetOperator{iDetOp}])
      if strcmp(Exp.DetOperator{iDetOp},'electron coherence') 
        LabelsDetectionOp{iDetOp} = 'signal (arb.u.)';
      else
        LabelsDetectionOp{iDetOp} = ['\langle' Exp.DetOperator{iDetOp} '\rangle'];
      end
    else
      % title of detection operators that were provided in matrix form
      title(['Signal detected with detection operator no.' num2str(iDetOp)])
      LabelsDetectionOp{iDetOp} = 'expectation value';
    end
    hold on
  end
end

% Get labels for indirect dimensions and set up axes
if isfield(Exp,'nPoints')
  LabelsIndirectDims = cell(1,length(Exp.nPoints));
  AxesIndirectDims = cell(1,length(Exp.nPoints));
  for iDim = 1 : length(Exp.nPoints)
    CurrentDim = ['Dim' num2str(iDim)];
    if length(Exp.(CurrentDim){1,2}) == 1
      % axis and its label in case of linear increments
      AxesIndirectDims{iDim} = Exp.(CurrentDim){1,2}*(0:Exp.nPoints(iDim)-1);
      LabelsIndirectDims{iDim} = [num2str(CurrentDim) ' / \Delta' Exp.(CurrentDim){1,1}];
    else 
      % axis and its label in case of userdefined increments
      AxesIndirectDims{iDim} = 1:Exp.nPoints(iDim);
      LabelsIndirectDims{iDim} = [num2str(CurrentDim) ' / Data Points ' Exp.(CurrentDim){1,1}];
    end
  end
end

TransientLabel = 'time (\mus)';

% plotting single point detection
if Opt.SinglePointDetection
  logmsg(1,'  single point detection:');
  if nDataPoints == 1
    % if only a single datapoint was acquired, there is no point in
    % plotting it, instead the output is displayed in the console
    for iDetOp = 1 : nDetOps
      if ischar(Exp.DetOperator{iDetOp})
        string  = ['  Expectation value of ' Exp.DetOperator{iDetOp} ':   '];
      else
        string = ['  Expectation value of operator no.' num2str(iDetOp) ':   '];
      end
      disp([string num2str(Signal(1,iDetOp))]);
    end
  else
    if length(Exp.nPoints) > 2
      % more than two dimensions can not be displayed in a general
      % way
      logmsg(1,'  more than two indirect dimensions - stopping');
      disp('Unable to display more than two indirect dimensions.')
    elseif length(Exp.nPoints) == 1
      % one dimensional case
      logmsg(1,'  creating plot(s) for one indirect dimension');
      for iDetOp = 1 : nDetOps
        figure
        plot(AxesIndirectDims{1},real(squeeze(Signal(:,1,iDetOp))));
        xlabel(LabelsIndirectDims{1});
        ylabel(LabelsDetectionOp{iDetOp})
      end
    elseif length(Exp.nPoints) == 2
      % two dimensional case
      logmsg(1,'  creating plot(s) for two indirect dimensions');
      for iDetOp = 1 : nDetOps
        figure(iDetOp)
        surf(AxesIndirectDims{1},AxesIndirectDims{2},real(squeeze(Signal(:,:,1,iDetOp)')));
        shading flat
        xlabel(LabelsIndirectDims{1});
        ylabel(LabelsIndirectDims{2})
        view([41.2000 61.2])
      end
    end
  end
  
else
% Transient plotting

  if isfield(Exp,'nPoints') && length(Exp.nPoints) > 1
    % more than 1 indirect dimension in combination with transient
    % detection can not be displayed
    logmsg(1,'  more than one indirect dimension in combination with transient - can not plot this - stopping');
    disp('Unable to display more than two indirect dimensions in combination with transient detection.')
  else
    % plotting transients
    logmsg(1,'  plotting transient(s)');
    if nDataPoints == 1
      % plotting a single acquisition point
      for iDetOp = 1 : nDetOps
        figure;
        plot(TimeAxis,real(Signal(:,iDetOp)));
        xlabel(TransientLabel)
        ylabel(LabelsDetectionOp{iDetOp})
      end
    else
      if ~iscell(Signal)
        % plotting if signals are organized in a numeric array
        SignalSize = size(Signal);
        if length(Exp.DetOperator) == 1
          SignalSize(end+1) = 1;
        end
        % reshape such to three dimensions
        Signal = reshape(Signal,[nDataPoints SignalSize(end-1) SignalSize(end)]);
      end
      
      for iDetOp = 1 : nDetOps
        figure;
        if iscell(Signal)
          for iDataPoint = 1 : nDataPoints
            % plotting the cell arrays
            SecondAxis = ones(1,length(TimeAxis{iDataPoint}))*AxesIndirectDims{1}(iDataPoint);
            plot3(TimeAxis{iDataPoint},SecondAxis,real(squeeze(Signal{iDataPoint}(:,iDetOp))));
          end
        elseif size(TimeAxis,1) == 1
          % plotting the reshaped numeric arrays
          surf(TimeAxis,AxesIndirectDims{1},real(squeeze(Signal(:,:,iDetOp))));
          shading flat
        else
          % plotting for when detected events have the same length,
          % but not identical timeaxes
          for iDataPoint = 1 : nDataPoints
            % plotting the cell arrays
            SecondAxis = ones(1,size(TimeAxis,2))*AxesIndirectDims{1}(iDataPoint);
            plot3(TimeAxis(iDataPoint,:),SecondAxis,real(Signal(iDataPoint,:,iDetOp)));
          end
        end
        
        view([41.2000 61.2])
        xlabel(TransientLabel)
        ylabel(LabelsIndirectDims{1})
        zlabel(LabelsDetectionOp{iDetOp})
        
      end
    end
  end
end