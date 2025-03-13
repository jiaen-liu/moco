function [c,r0]=cubic_intp(r,nr,a)
    if nargin<3
        a=-0.5;
    end
    si=size(r);
    ndim=si(1);
    n=si(2);
    r0=zeros(4,ndim,n);
    c=zeros(4,ndim,n);
    % interpolation locations
    for i=1:ndim
        r0(2,i,:)=floor(r(i,:));
        r0(1,i,:)=r0(2,i,:)-1;
        r0(3,i,:)=r0(2,i,:)+1;
        r0(4,i,:)=r0(2,i,:)+2;
        % && extend boundary
        % The kernel is in 0~nr-1
        for j=1:size(r0,1)
            m=r0(j,i,:)<0;
            r0(j,i,m)=0;
            m=r0(j,i,:)>nr(i)-1;
            r0(j,i,m)=nr(i)-1;
        end
    end
    % interpolation values
    for i=1:ndim
        c(:,i,:)=cubic_kernel(mod(r(i,:),1),a);
    end
end

function y=cubic_kernel(x,a)
% x: normalized to dx=1
    if nargin<2
        a=-0.5;
    end
    x=x(:);
    nx=length(x);
    y=zeros(4,nx);
    y(1,:)=cubic_spline(-1-x,a);
    y(2,:)=cubic_spline(-x,a);
    y(3,:)=cubic_spline(1-x,a);
    y(4,:)=cubic_spline(2-x,a);
end

function y=cubic_spline(x,a)
    y=zeros(size(x));
    m=abs(x) <= 1;
    y(m)=(a+2.0)*abs(x(m)).^3.0-(a+3.0)*x(m).^2+1.0;
    m=abs(x)<=2 & abs(x)>1;
    y(m)=a*(abs(x(m)).^3-5.0*x(m).^2+8.0*abs(x(m))-4.0);
end
