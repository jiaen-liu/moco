function [irow,icol,val]=interp_kernel(coord0,coord1,varargin)
    si0=size(coord0);
    si1=size(coord1);

    method='linear';

    nvarargin=numel(varargin);
    if nvarargin>0
        method=varargin{1};
    end
    if numel(si0)==3 && si0(3)==2
        ndim=2;
    elseif numel(si0)==4 && si0(4)==3
        ndim=3;
    else
        error('The input coordinate should have two or three dimensions!');
    end
    nx1=si1(1);
    ny1=si1(2);
    nz1=1;
    nx0=si0(1);
    ny0=si0(2);
    nz0=1;
    if numel(si1)==4
        nz1=si1(3);
    end
    if numel(si0)==4
        nz0=si0(3);
    end
    n1=nx1*ny1*nz1;
    % normalize coordinate 1 according to coordinate 0
    [x1,x2,x3]=norm_coord(coord0,coord1);
    switch method
      case 'linear'
        nker=2^ndim;
      case 'cubic'
        nker=4^ndim;
    end
    % row index
    irow=repmat([1:n1],[nker,1]);
    irow=irow(:);
    % column index
    icol=zeros(nker,n1);
    % value of each elements
    val=zeros(nker,n1);
    %
    if ndim==2
        r=[x1(:),x2(:)]';
        nr=[nx0,ny0];
    elseif ndim==3
        r=[x1(:),x2(:),x3(:)]';
        nr=[nx0,ny0,nz0];
    end
    % r0 is zero-based
    switch method
      case 'cubic'
        [c,r0]=cubic_intp(r,nr,-0.5);
      case 'linear'
        [c,r0]=linear_intp(r,nr);
    end
    sic=size(c);
    ks=sic(1);
    switch ndim
      case 2
        for i=1:ks
            for j=1:ks
                icol(i+(j-1)*ks,:)=r0(i,1,:)+1+(r0(j,2,:))*nx0;
                val(i+(j-1)*ks,:)=c(i,1,:).*c(j,1,:);
            end
        end
      case 3
        for i=1:ks
            for j=1:ks
                for k=1:ks
                    icol(i+(j-1)*ks+(k-1)*ks^2,:)=...
                        r0(i,1,:)+1+r0(j,2,:)*nx0+r0(k,3,:)*(nx0*ny0);
                    val(i+(j-1)*ks+(k-1)*ks^2,:)=...
                        c(i,1,:).*c(j,2,:).*c(k,3,:);
                end
            end
        end
    end
    % irow and icol are one-based
    irow=irow(:);
    icol=icol(:);
    val=val(:);
end