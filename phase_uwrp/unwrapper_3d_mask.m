function pha_uw=unwrapper_3d_mask(pha,mask)
    pha(isnan(pha))=0;
    pha=single(pha);
    if nargin < 2
        mask=ones(size(pha));
    end
    mask_cp=zeros(size(mask),'uint8');
    mask_cp(mask~=0)=255;
    if length(size(pha))==2
        [nx,ny]=size(pha);
        pha=reshape(pha,[nx,ny,1]);
        mask_cp=reshape(mask_cp,[nx,ny,1]);
    end
    pha_uw=double(unwrapper_3d_mask_matlab(pha,mask_cp));
end