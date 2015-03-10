function [avg, asympAnalysis, bwX, bwY] = getTwoDistroInfFunAvgs(X, Y, ...
  infFunX, infFunY, asympVarFun, params)
% This function will sum the influence functions over the partitions

  numPartitions = params.numPartitions;
  infFunXPartTerms = zeros(numPartitions, 1);
  infFunYPartTerms = zeros(numPartitions, 1);
  asympVarXPartTerms = zeros(numPartitions, 1);
  asympVarYPartTerms = zeros(numPartitions, 1);
  partWeightsX = zeros(numPartitions, 1);
  partWeightsY = zeros(numPartitions, 1);
  n = size(X, 1);
  m = size(Y, 1);

  for k = 1:numPartitions

    [Xden, Xest] = getDenEstSamples(X, numPartitions, k);
    [Yden, Yest] = getDenEstSamples(Y, numPartitions, k);

      % First determine the bandwidths for density estimation
      if k == 1
        % For X data
        if ~isfield(params, 'bandwidthX') | isempty(params.bandwidthX)
          bwX = kdePickBW(Xden, params.smoothness, params);
        else
          bwX = params.bandwidthX(Xden);
        end
        % For Y data
        if ~isfield(params, 'bandwidthY') | isempty(params.bandwidthY)
          bwY = kdePickBW(Yden, params.smoothness, params);
        else
          bwY = params.bandwidthY(Yden);
        end
      end

    % Obtain the KDE at
%     size(Xden), size(Xest), size(Yden), size(Yest),
    densEstX = kdeGivenBW(Xden, bwX, params.smoothness, params);
    densEstY = kdeGivenBW(Yden, bwY, params.smoothness, params);
    densXatX = densEstX(Xest);
    densXatY = densEstX(Yest);
    densYatX = densEstY(Xest);
    densYatY = densEstY(Yest);

    % Now obtain the sum of influence functions
    infFunXPartTerms(k) = sum( infFunX(densXatX, densYatX) );
    infFunYPartTerms(k) = sum( infFunY(densXatY, densYatY) );

    % If doing asymptotic analysis
    if params.doAsympAnalysis
      partWeightsX(k) = size(Xest, 1);
      partWeightsY(k) = size(Yest, 1);
      [aX, aY] = asympVarFun(densXatX, densXatY, densYatX, densYatY);
      asympVarXPartTerms(k) = aX;
      asympVarYPartTerms(k) = aY;
    end

  end

  % Now return the average as the influence function
  avg = sum(infFunXPartTerms)/n + sum(infFunYPartTerms)/m;

  if params.doAsympAnalysis
    asympVarX = (partWeightsX' * asympVarXPartTerms)/n;
    asympVarY = (partWeightsY' * asympVarYPartTerms)/m;
    asympVar = (n+m) * (asympVarX/n + asympVarY/m);
    asympStd = sqrt(asympVar);
    asympAnalysis.asympVar = asympVar;
    asympAnalysis.asympStd = asympStd;

    % Now construct the confidence interval
    if isfield(params, 'alpha')
      w = norminv(1-params.alpha/2);
      % If avg is not the estimate then you need to correct this bias later
      asympAnalysis.confInterval(1) = avg - w * asympStd/sqrt(n+m);
      asympAnalysis.confInterval(2) = avg + w * asympStd/sqrt(n+m);
    end

  else
    asympAnalysis = [];
  end
  
end
