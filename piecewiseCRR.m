function price = piecewiseCRR(S0, K, r, T, sigma, divs, divt, N)
% Binomial Tree Model of CRR with discrete dividends
dt = T / N;
u = exp(sigma * sqrt(dt));
d = 1 / u;
p = (exp(r * dt) - d) / (u - d);
discount = exp(- r * dt);
ddiv = ceil(divt * N / T) + 1;
divl = zeros(1,N+1);
count = 0;
for i = 1:(N+1)
    if (i <= ddiv(1))
        svals = zeros(1,i);
        divl(i) = i;
        for j = 1:i
            svals(j) = S0 * u^(i - 1) * (d/u)^(j - 1);
        end
    else
        len = i - ddiv(count) + 1;
        num = divl(ddiv(count));
        ancester = ddiv(count);
        ii = len * num;
        divl(i) = ii;
        svals = zeros(1, ii);
        for j = 1:ii
            svals(j) = SVals{ancester}(floor((j-1)/len)+1) * u^(len-1) * (d/u)^mod(j-1,len);
        end   
    end
    if (count < sum(divs>=0))
        if (i == ddiv(count+1))
            svals = max(svals - divs(count+1),0);
            count = count + 1;
        end
    end
    SVals{i} = svals;
end
PVals{N + 1} = max(K - SVals{N + 1}, 0);
for i = N:-1:1
    pvals = zeros(1, divl(i));
    if (count && i >= ddiv(count))
        len = i - ddiv(count) + 1;
        if (i == ddiv(count))
            count = count - 1;
        end
    end
    for j = 1:divl(i)
        if (i >= ddiv(1) && mod(j-1,len) == 0)
            up = floor((j-1)/len) * (len + 1) + 1;
        end
        if (i < ddiv(1))
            up = j;
        end
        Early = discount * (p * PVals{i + 1}(up) + (1 - p) * PVals{i + 1}(up + 1));
        up = up + 1;
        Internal = max(K - SVals{i}(j), 0);
        pvals(j) = max(Early, Internal);
    end
    PVals{i} = pvals;
end
price = PVals{1};
% price = BinomialTreeCRRDiscrete(50,50,0.05,5/12,0.4,[2],[2/12],1000)
end