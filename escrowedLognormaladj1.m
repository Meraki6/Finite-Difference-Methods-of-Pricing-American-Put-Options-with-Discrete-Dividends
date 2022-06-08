function price = escrowedLognormaladj1(S0, K, r, T, sigma, divs, divt, N)
% Binomial Tree Model of lognormal with discrete dividends
dt = T / N;
S0d = S0 - sum(divs .* exp(- r * divt));
sigma = sigma * S0 / S0d;
u = exp((r - 0.5 * sigma ^ 2) * dt + sigma * sqrt(dt));
d = exp((r - 0.5 * sigma ^ 2) * dt - sigma * sqrt(dt));
p = 0.5;
discount = exp(- r * dt);
ddiv = ceil(divt * N / T) + 1;
for i = 1:(N+1)
    svals = zeros(1, i);
    for j = 1:i
        svals(j) = S0d * u^(i - 1) * (d/u)^(j - 1) + sum(divs .* exp(-r*dt*(ddiv-i)) .* (i<ddiv));
    end
    SVals{i} = svals;
end
PVals{N + 1} = max(K - SVals{N + 1}, 0);
for i = N:-1:1
    pvals = zeros(1, i);
    for j = 1:i
        Early = discount * (p * PVals{i + 1}(j) + (1 - p) * PVals{i + 1}(j + 1));
        Internal = max(K - SVals{i}(j), 0);
        pvals(j) = max(Early, Internal);
    end
    PVals{i} = pvals;
end
price = PVals{1};
% price = BinomialTreeLognormalDiscrete(50,50,0.05,5/12,0.4,[2],[2/12],1000)
end

