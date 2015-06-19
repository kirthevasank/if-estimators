function [avg, asympAnalysis, bw] = getInfFunAvgs(X, infFun, asympVarFun,params)
% This function will sum the influence functions over the partitions.

  n = size(X, 1);
  numPartitions = params.numPartitions;
  numAvgPartitions = params.numAvgPartitions;

  % Now determine the number of partitions to average over.
  infFunPartTerms = zeros(numAvgPartitions, 1);
  asympVarPartTerms = zeros(numAvgPartitions, 1);
  partWeights = zeros(numAvgPartitions, 1);

  for k = 1:numPartitions

    [Xden, Xest] = getDenEstSamples(X, numPartitions, k);

      % First determine the bandwidth for density estimation
      if k == 1 
        if ~isfield(params, 'bandwidth') | isempty(params.bandwidth)
          bw = kdePickBW(Xden, params.smoothness, params);
        else
          bw = params.bandwidth(Xden);
        end
      end

    % Obtain the KDE
    densEst = kdeGivenBW(Xden, bw, params.smoothness, params);
    densX = densEst(Xest);

    % Now obtain the sum of influence function values for Xest
    infFunPartTerms(k) = sum( infFun(densX) );
    partWeights(k) = size(Xest, 1);

    % If doing asymptotic analysis
    if params.doAsympAnalysis
      asympVarPartTerms(k) = asympVarFun(densX);
    end

  end

  % Now return the average
  avg = sum(infFunPartTerms)/sum(partWeights);

  if params.doAsympAnalysis
    asympVar = (partWeights' * asympVarPartTerms)/n;
    asympStd = sqrt(asympVar);
    asympAnalysis.asympVar = asympVar;
    asympAnalysis.asympStd = asympStd;
    if asympVar < 0
      fprintf(['The estimated asymptotic variance is negative. This is', ...
        ' probably because the asymptotic distribution is degenerate and ', ...
         'non-gaussian.\n']);
    end

    % Now construct a confidence interval
    if isfield(params, 'alpha')
      w1 = norminv(1 - params.alpha/2);
      % If avg is not the estimate then you need to correct the bias later !
      asympAnalysis.confInterval(1) = avg - w1 * asympStd/sqrt(n);
      asympAnalysis.confInterval(2) = avg + w1 * asympStd/sqrt(n);
    end

  else
    asympAnalysis = [];
  end

end

