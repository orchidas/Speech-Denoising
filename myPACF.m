function [ pacf ] = myPACF(x, p)
%calculates pacf from time series signal x and model order p using
%Yule-Walker equations

pacf = zeros(p,1);
% get autocorrelation coefficients
[autocorr, lags] = xcorr(x , 'coeff');

for k = 1:p
    R =zeros(k,k);
    r = autocorr(find(lags == 1) : find(lags == k))';
    for i = 1:k
     for j = 1:k
         R(i,j) = autocorr(lags == i-j);
     end
    end
    Phi = R\r; %inverse(R) * r
    pacf(k) = Phi(k);
end


end

