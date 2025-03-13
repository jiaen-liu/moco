function y=fit_pha0(freq,t,x)
t=t(:);
x=x(:);
m=abs(x);
p=angle(x);
dp=freq*t*2*pi-p;
a1=(m.^2).'*cos(dp);
a2=-(m.^2).'*sin(dp);
ma=(a1^2+a2^2)^0.5;
sin_y=a1/ma;
cos_y=a2/ma;
if sin_y~=0
    beta=acos(cos_y)*sign(sin_y);
elseif cos_y>=0;
    beta=0;
else
    beta=-pi;
end
y=pi/2-beta;
end