function h = fgaussian3(siz,sig)
    if numel(siz)==1
        siz=siz*ones(3,1);
    end
    if numel(sig)==1
        sig=sig*ones(3,1);
    end
    siz   = floor((siz-1)/2);
    [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
    h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2));
    h = h/sum(h(:));       
end