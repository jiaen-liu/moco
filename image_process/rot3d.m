function m=rot3d(thetaz, thetay, thetax)
% $$$   rotation depends on the order
% $$$   for the function, it starts with x, y and ends on z
% $$$   applies to right hand coordinate frame and right hand rotation or
% $$$   left hand coordi and left hand rotation
    if nargin == 1
        thetax=thetaz(1);
        thetay=thetaz(2);
        thetaz=thetaz(3);
    end
    thetaz = double(thetaz);
    thetay = double(thetay);
    thetax = double(thetax);
    m=[[cos(thetaz)*cos(thetay), sin(thetaz)*cos(thetay), -sin(thetay)].', ...
       [-sin(thetaz)*cos(thetax)+cos(thetaz)*sin(thetay)*sin(thetax), cos(thetaz)*cos(thetax)+sin(thetaz)*sin(thetay)*sin(thetax), cos(thetay)*sin(thetax)].', ...
       [sin(thetaz)*sin(thetax)+cos(thetaz)*sin(thetay)*cos(thetax), -cos(thetaz)*sin(thetax)+sin(thetaz)*sin(thetay)*cos(thetax), cos(thetay)*cos(thetax)].'];
end