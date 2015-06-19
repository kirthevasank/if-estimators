% This is a demo showing the simplest use case.
% Also see demo1, ..., demo4 for different functionals and obtaining asymptotic
% confidence sets.

X = rand(1000, 2);
Y = rand(1000, 2);

% Estimate the Shannon Entropy
% The empty matrices are for default settings.
estim = shannonEntropy(X, [], []); 
fprintf('Shannon Entropy: Estimate: %0.4f, Truth: 0\n', estim);

% Estimate the Hellinger Divergence
estim = hellingerDivergence(X, Y, [], []); 
fprintf('Hellinger Divergence: Estimate: %0.4f, Truth: 0\n', estim);

