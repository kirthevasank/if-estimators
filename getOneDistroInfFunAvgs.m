function [avg, asympAnalysis, bw] = getOneDistroInfFunAvgs(...
  X, infFun, asympVarFun, params)
% This function will sum the influence functions over the partitions.

  numPartitions = params.numPartitions;
  infFunPartTerms = zeros(numPartitions, 1);
  asympVarPartTerms = zeros(numPartitions, 1);
  partWeights = zeros(numPartitions, 1);
  n = size(X, 1);

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

    % If doing asymptotic analysis
    if params.doAsympAnalysis
      asympVarPartTerms(k) = asympVarFun(densX);
      partWeights(k) = size(Xest, 1);
    end

  end

  % Now return the average
  avg = sum(infFunPartTerms)/n;

  if params.doAsympAnalysis
    asympVar = (partWeights' * asympVarPartTerms)/n;
    asympStd = sqrt(asympVar);
    asympAnalysis.asympVar = asympVar;
    asympAnalysis.asympStd = asympStd;

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

