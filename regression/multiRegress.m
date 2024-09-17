function [c,y,resi]=multiRegress(d,X,mask)
    n=size(X,1);
    if nargin<3
        mask=true(n,1);
    end
    c=X(find(mask),:)\d(find(mask),:);
    if nargout>1
        y=X*c;
    end
    if nargout>2
        resi=d-y;
    end
end