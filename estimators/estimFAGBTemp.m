function tempEstim= estimFAGBTemp(densXatX, densXatY, densYatX, densYatY, alpha)
% An estimate for \int f^alpha g^beta to be used for the asymptotic variance in
% several functionals (e.g: fAlphaGBeta, KL)
  beta = 1-alpha;
  tempEstim = alpha * mean( (densYatX ./ densXatX).^beta ) + ...
              beta * mean( (densXatY ./ densYatY).^alpha );
end

