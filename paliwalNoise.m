function [ R ] = paliwalNoise(x, p)
%Paliwal's method of measurement noise calculation

[c,lags]=xcorr(x);
si=0;
so=0;
A=lpc(x,p);

%keeping last len/2 samples of autocorrelation as it is an even function of
%lags
c = c(find(lags == 0):length(c));
lags = lags(find(lags == 0):length(lags));

for i = 1:p 
    for k=1:p
        si=si+ A(k)*c(lags == abs(i-k));
    end
    so = so + (A(i)*(c(lags == i)+si));
end

R = so/sum(A.^2)


% s = 0;
% for i = 1:p
%     s = s + A(i)*c(lags == i);
% end
% 
% R = s + c(lags == 0)

end

