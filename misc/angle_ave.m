function y=angle_ave(a)
% radian, -pi to pi
a=a(:);
if max(a)-min(a)<pi
    amin=min(a);
    a=a-amin;
    y=mean(a)+amin;
else
    a=wrapTo2Pi(a);
    amin=min(a);
    a=a-amin;
    y=mean(a)+amin;
    y=wrapToPi(y);
end
end