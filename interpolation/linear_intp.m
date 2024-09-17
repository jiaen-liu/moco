function [c,r0]=linear_intp(r,nr)
    si=size(r);
    ndim=si(1);
    n=si(2);
    r0=zeros(2,ndim,n);
    c=zeros(2,ndim,n);
    % interpolation locations
    for i=1:ndim
        r0(1,i,:)=floor(r(i,:));
        r0(2,i,:)=r0(1,i,:)+1;
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
        c(:,i,:)=linear_kernel(mod(r(i,:),1));
    end
end

function y=linear_kernel(x)
    x=x(:);
    n=length(x);
    y=zeros(2,n);
    y(1,:)=1-x;
    y(2,:)=x;
end


