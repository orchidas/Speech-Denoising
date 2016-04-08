function [Y] = adjustDimensions(X,p)
%Adjust the dimensions of X to pxp
m = size(X,1);
Y = zeros(p,p);
if(p > m)
    newRows = zeros(p-m,m);
    newCols = zeros(p,p-m);
    temp = [X;newRows];
    Y = [temp newCols];
else
    if(p == m)
        Y = X;
    else    
        Y(:,:) = X(1:p,1:p);
    end
end
end

