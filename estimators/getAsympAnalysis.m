function asympAnalysis = getAsympAnalysis(estim, asympVar, alpha, n)
% Returns a structure containing some statistics for constructing asymptotic
% confidence intervals.
  asympAnalysis.asympVar = asympVar;
  asympAnalysis.asympStd = sqrt(asympVar);
  width = norminv(1 - alpha/2) * asympAnalysis.asympStd/sqrt(n);
  asympAnalysis.confInterval(1) = estim - width;
  asympAnalysis.confInterval(2) = estim + width;
end

