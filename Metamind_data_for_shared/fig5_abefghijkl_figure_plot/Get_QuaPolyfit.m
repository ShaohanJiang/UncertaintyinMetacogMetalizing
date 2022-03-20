function P = Get_2rdPolyfit(tmpa)

[nsub, nx]= size(tmpa);
P = nan(nsub,3);
for jj = 1:nsub
    x = 1:nx;
    y = tmpa(jj, :);
    x = zscore(x);
    y = zscore(y);
    n = 2;
    p = polyfit(x, y, n);
    P(jj, :) = p;
    
end