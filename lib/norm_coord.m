function [x1,x2,x3]=norm_coord(coord0,coord1)
    si0=size(coord0);
    x1_0=si0(1);
    x2_0=si0(2);
    if numel(si0) == 3
        x3_0 = 1;
        coord0 = reshape(coord0, [x1_0, x2_0, x3_0, 3]);
    else
        x3_0 = si0(3);
    end;

    si1=size(coord1);
    x1_1=si1(1);
    x2_1=si1(2);
    if numel(si1) == 3
        x3_1 = 1;
        coord1 = reshape(coord1, [x1_1, x2_1, x3_1, 3]);
    else
        x3_1 = si1(3);
    end;
    dx1_0=reshape(coord0(2,1,1,:)-coord0(1,1,1,:),[3,1]);
    dx2_0=reshape(coord0(1,2,1,:)-coord0(1,1,1,:),[3,1]);
    

    dx1_0_squ=sum(abs(dx1_0).^2);
    dx2_0_squ=sum(abs(dx2_0).^2);

    
    n=x1_1*x2_1*x3_1;
    x1=zeros(n,1);
    x2=zeros(n,1);
    x3=zeros(n,1);

    coord1_rfm=reshape(coord1,[n,3]);

    x1=(coord1_rfm*dx1_0-reshape(coord0(1,1,1,:),[1,3])*dx1_0)/dx1_0_squ;
    x2=(coord1_rfm*dx2_0-reshape(coord0(1,1,1,:),[1,3])*dx2_0)/dx2_0_squ;

    if x3_0 > 0
        dx3_0=reshape(coord0(1,1,2,:)-coord0(1,1,1,:),[3,1]);
        dx3_0_squ=sum(abs(dx3_0).^2);
        x3=(coord1_rfm*dx3_0-reshape(coord0(1,1,1,:),[1,3])*dx3_0)/dx3_0_squ;
    end
end