% This is a demo showing the simplest use case for all functionals.
% Also see demo1, ..., demo4 for different functionals and obtaining asymptotic
% confidence sets.
% The empty matrices are for default settings.

X = rand(2000, 1);
Y = rand(2000, 1);
Z = rand(2000, 1);
W = rand(2000, 1);
fp.alpha = 0.75; % for renyi/ tsallis divergences

% Entropies
fprintf('Entropies\n-----------------------------------------------------\n');

estim = shannonEntropy(X, [], []); 
fprintf('Shannon Entropy: Estimate: %0.4f, Truth: 0\n', estim);


% Divergences
fprintf('\nDivergences\n--------------------------------------------------\n');

estim = hellingerDivergence(X, Y, [], []); 
fprintf('Hellinger Divergence: Estimate: %0.4f, Truth: 0\n', estim);

estim = klDivergence(X, Y, [], []); 
fprintf('KL Divergence: Estimate: %0.4f, Truth: 0\n', estim);

estim = chiSqDivergence(X, Y, [], []); 
fprintf('Chi-Squared Divergence: Estimate: %0.4f, Truth: 0\n', estim);

estim = tsallisDivergence(X, Y, fp, []); 
fprintf('Tsallis Divergence: Estimate: %0.4f, Truth: 0\n', estim);

estim = renyiDivergence(X, Y, fp, []); 
fprintf('Renyi Divergence: Estimate: %0.4f, Truth: 0\n', estim);


% Mutual Informations
fprintf('\nMutual Informations\n------------------------------------------------\n');

estim = shannonMI(X, Y, [], []); 
fprintf('Shannon MI: Estimate: %0.4f, Truth: 0\n', estim);


% Conditional Quantities 
fprintf('\nConditional Quantities\n---------------------------------------------\n');

estim = condShannonEntropy(X, Z, [], []); 
fprintf('Cond. Shannon Entropy: Estimate: %0.4f, Truth: 0\n', estim);

estim = condKLDivergence(X, Y, Z, W, [], []); 
fprintf('Cond. KL Divergence: Estimate: %0.4f, Truth: 0\n', estim);

estim = condTsallisDivergence(X, Y, Z, W, fp, []); 
fprintf('Cond. Tsallis Divergence: Estimate: %0.4f, Truth: 0\n', estim);

estim = condShannonMI(X, Y, Z, [], []); 
fprintf('Cond. Shannon Divergence: Estimate: %0.4f, Truth: 0\n', estim);

