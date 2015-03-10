function utX
% Unit tests for all one distro functionals

  close all;
  clear all;

  functionals = {'shannonEntropy'};
  tests = {'1D-Uniform', '1D-Conv', '2D-Conv', '2D-Gaussian'};
  functionalParams = struct;
  params = struct;
  params.alpha = 0.05;
  params.doAsympAnalysis = true;
  params.kdePickMethod = 'silverman';
  params.kdePickMethod = 'cv';

  % Test 1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for testIdx = 1:numel(tests)

    % First generate the data
    if testIdx == 1
      params.doBoundaryCorrection = true;
      N = 1000;
      X = rand(N, 1);
      trueVals = [0];

    elseif testIdx == 2
      params.doBoundaryCorrection = true;
      N = 1000;
      gamma = 10;
      Z = rand(N, 1+gamma); B = double(rand(N, 1) < 0.5);
      X = B.* Z(:,1) + (1-B).*max(Z(:,2:end), [], 2);
      trueDensity = @(t) 0.5 + 0.5*gamma* t.^(gamma-1);
      params.estLowerBound = 0.4;

      % Compute the true Entropy
      entropyFunc = @(t) trueDensity(t) .* log( trueDensity(t) );
      t = linspace(0,1,1000); trueEntropy = -mean(entropyFunc(t));
      
      % True vals
      trueVals(1) = trueEntropy;

    elseif testIdx == 3
      params.doBoundaryCorrection = true;
      N = 1000;
      gamma = 10;
      Z = rand(N, 1+gamma); B = double(rand(N, 1) < 0.5);
      Y = B.* Z(:,1) + (1-B).*max(Z(:,2:end), [], 2);
      X = [Y rand(N, 1)];
      params.estLowerBound = 0.4;
      % Since the 2nd dimension is the uniform density, the trueVals will be the
      % same as for testIdx = 2

    elseif testIdx == 4
      params.doBoundaryCorrection = false;
      d = 5;
      N = 10000;
      X = randn(N, d);
      trueVals(1) = d/2 *(1 + log(2*pi));
      params.estLowerBound = 0;

    end

    % Now do the tests
    doTests(X, functionals, trueVals, tests{testIdx}, params, functionalParams);

  end

end


function doTests(X, functionals, trueVals, test, params, functionalParams)

  fprintf('%s\n====================================================\n\n', test);

  for k = 1:numel(functionals)

    functional = functionals{k};
    trueVal = trueVals(k);
    fprintf('  %s: Truth: %0.5f\n', functional, trueVal);

    % First Do a Data split version
%     params.dataSplit = true;
    if strcmp(functional, 'shannonEntropy')
      [estDS, asympAnalysis] = shannonEntropy(X, functionalParams, params);
    end
    % Now compute the errors
    errDS = abs(trueVal - estDS);
    fprintf('    EstimDS : %.4f,  ErrDS : %.4f, CI: %s\n', ...
      estDS, errDS, mat2str(asympAnalysis.confInterval));

%     % Now do a non-Data split version
%     params.dataSplit = false;
%     [estNDS, asympAnylsis] = estimateOneDistroFunctionals...
%       (X, functional, functionalParams, params);
%     % Now compute the errors
%     errNDS = abs(trueVal - estNDS);
%     if trueVal ~= 0
%       errNDS = errNDS/trueVal;
%     end
%     fprintf('    EstimNDS : %.4f,  ErrDS : %.4f, CI: %s\n', ...
%       estNDS, errNDS, mat2str(asympAnylsis.confInterval));
% 
%     fprintf('\n');
  end

end

