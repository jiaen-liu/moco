function w = win_tukey (L, r, x)
    if nargin<2
        r = 1/2;
    end
    if nargin<3
        x=[0:L-1]'/(L-1);
    end
    if any(x<0 | x>1)
        error('*** The explanatory variable x should be in [0,1]! ***');
    end
    w = zeros(size(x));
    w(x >= r/2.0 & x < (1.0-r/2.0)) = 1.0;
    idx = find(x >= 0 & x < r/2.0);
    w(idx) = 0.5*(1.0+cos(2.0*pi/r*(x(idx)-r/2.0)));
    idx = find(x >= 1.0-r/2.0);
    w(idx) = 0.5*(1.0+cos(2.0*pi/r*(x(idx)-1.0+r/2.0)));
end