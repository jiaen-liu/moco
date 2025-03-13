function y=volerode(im,n)
    si=size(im);
    nx=si(1);
    ny=si(2);
    if numel(si)==3
        nz=si(3);
    else
        nz=1;
    end
    y=zeros(si);
    se=strel('disk',floor(n));
    for i=1:nz
        y(:,:,i)=imerode(im(:,:,i),se);
    end
end