function y=ref_mag(x,mask,frac)
    imt=sum(abs(x),4);
    [nx,ny,nz,nch]=size(x);
    if nnz(imt)==0
        y=zeros(nx,ny,nz);
        return;
    end
    % smooth the image
    ims=smoes(imt,frac,mask,1);
    y=imt./(max(ims, max(ims(:))/2.5));
end