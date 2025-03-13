function [sensit,mask]=interp_sense(coord0,coord1,sensit,mask)
    sensit=interp3_nmat(coord0,coord1,sensit,'linear','none');
    sensit(isnan(sensit))=0;
    if nargin>3
        mask=interp3_nmat(coord0,coord1,double(mask),'linear','none');
        mask(isnan(mask))=0;
        mask=logical(mask);
    else
        mask=[];
    end
end