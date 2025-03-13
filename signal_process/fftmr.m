function y=fftmr(data,direction,dim)

if numel(dim) > 1
    y=data;
    for i=1:numel(dim)
        y=fftmr(y,direction,dim(i));
    end
    return
end
si=size(data);
ndim=numel(si);
if dim > ndim
    y=data;
else
    if direction == 1 
        y=si(dim)*fftshift(ifft(ifftshift(data,dim),si(dim),dim),dim);
    elseif direction == -1
        y=fftshift(fft(ifftshift(data,dim),si(dim),dim)/si(dim),dim);
    end
end
end