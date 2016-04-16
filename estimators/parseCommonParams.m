function params = parseCommonParams(params, numDims, n)
% params is a struct with some parameters for the estimation. 
% (viz. whether to data split or not)
% numDims is the number of Dimensions

  if ~exist('params', 'var') | isempty(params)
    params = struct;
  end

  % Whether to do Asymptotic Analysis
  if ~isfield(params, 'doAsympAnalysis')
    params.doAsympAnalysis = false;
  end

  % The smoothness of the function for the KDE
  if ~isfield(params, 'smoothness')
      params.smoothness = 'gauss';
  end

  % Number of partitions to split the data into
  if ~isfield(params, 'numPartitions')
    params.numPartitions = 1; % by default, do not partition the data.
  end
  if isstr(params.numPartitions) & strcmp(params.numPartitions, 'loo')
    params.numPartitions = n;
  end

  % Number of partitions to average over
  if ~isfield(params, 'numAvgPartitions')
    if isfield(params, 'averageAll') & (~params.averageAll)
      params.numAvgPartitions = 1;
    else
      params.numAvgPartitions = params.numPartitions;
    end
  end

  % Some parameters for Kernel Density estimation
  if ~isfield(params, 'doBoundaryCorrection')
      params.doBoundaryCorrection = false;
  end

  if ~isfield(params, 'estLowerBound')
    params.estLowerBound = 1e-5;
  end

  if ~isfield(params, 'estUpperBound')
    params.estUpperBound = Inf;
  end

end

