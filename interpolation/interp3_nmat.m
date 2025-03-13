function y=interp3_nmat(coord0,coord1,v,varargin)
    si0=size(coord0);
    si1=size(coord1);
    siv=size(v);
    if numel(siv)<3
        siv=[siv,1,1];
    end
    method='linear';
    ex_method='linear';
    nvarargin=numel(varargin);
    if nvarargin>0
        method=varargin{1};
    end
    if nvarargin>1
        ex_method=varargin{2};
    end
    if numel(si0) ~= 4 || numel(si1) ~= 4 || ~isequal(si0(1:3), siv(1:3))
        error('The input coordinate should have three dimensions!');
    end
    nv=numel(v)/prod(siv(1:3));
    v=reshape(v,[siv(1:3),nv]);
    nx1=si0(1);
    nx2=si0(2);
    nx3=si0(3);
    [x1,x2,x3]=norm_coord(coord0,coord1);
    [x10,x20,x30]=ndgrid([0:nx1-1],[0:nx2-1],[0:nx3-1]);
    y=zeros(si1(1),si1(2),si1(3),nv);
    for i=1:nv
        F=griddedInterpolant(x10,x20,x30,v(:,:,:,i),method,ex_method);
        y(:,:,:,i)=reshape(F(x1,x2,x3),[si1(1),si1(2),si1(3)]);
    end
    y=reshape(y,[si1(1),si1(2),si1(3),siv(4:end)]);
end