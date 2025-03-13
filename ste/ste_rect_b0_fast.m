function db0=ste_rect_b0_fast(db0,par)
% $$$     db0_uc=db0(:,:,:,:,~fcrptd);
% $$$     db0_uc_ave=mean(db0_uc,4);
% $$$     db0_uc_dif=mean(db0_uc-db0_uc_ave,5);
% $$$     db0=db0-db0_uc_dif;
    if par.rect_cyc==0
        db0_uc=db0;
        db0_uc_ave=mean(db0_uc,4);
        db0_uc_dif=median(db0_uc-db0_uc_ave,5);
        db0=db0-db0_uc_dif;
    elseif par.rect_cyc==1
        si=size(db0);
        cyc=si(4);
        db0=combine_dim(db0,[1,2,3]);
        db0=combine_dim(db0,[2,3]);
        mask=false(si(4)*si(5),1);
        mask(si(4)+1:si(4)*(si(5)-1))=true;
        db0=steRegress(db0,cyc,cyc/2,mask);
        db0=reshape(db0,si);
    end
end