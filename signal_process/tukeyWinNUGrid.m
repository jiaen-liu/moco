function w=tukeyWinNUGrid(n, r, x)
  w = zeros(n, 1);
  if nargin<3
      x=[0:n-1].'/(n-1);
  end
  % asymmetric
  w(find(x >= r/2.0 & x < 1.0-r/2.0)) = 1.0;
  idx = find(x >= 0 & x < r/2.0);
  w(idx) = 0.5*(1.0+cos(2.0*pi/r*(x(idx)-r/2.0)));
  idx = find(x >= 1.0-r/2.0);
  w(idx) = 0.5*(1.0+cos(2.0*pi/r*(x(idx)-1.0+r/2.0)));
end
