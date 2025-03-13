function ste_info=ste_discard(ste_info,sorted_ste)
% discard data in motion-corrupted full-fov volumns
    nshot=sorted_ste.nshot;
    idx_shot_cenff=sorted_ste.idx_shot_cenff;
    idx_shot_cenpf=sorted_ste.idx_shot_cenpf;
    fcrptd=ste_info.fcrptd;
    mask_gre=logical(interp1(idx_shot_cenff,...
                             double(~fcrptd(1:numel(idx_shot_cenff))),...
                             [1:nshot],...
                             'nearest','extrap'));
    nk_shot=sorted_ste.para.nk_shot;
    for i=1:ste_info.nc
        idx_gre_cl=ste_info.idx_gre==i;
        mask_gre_cl=mask_gre(idx_gre_cl);
        ste_info.md{i}.nk=total(mask_gre_cl)*nk_shot;
        ktmp=ste_info.md{i}.k;
        ktmp=reshape(ktmp,[2,nk_shot,numel(ktmp)/2/nk_shot]);
        ktmp=ktmp(:,:,find(mask_gre_cl));
        ste_info.md{i}.k=combine_dim(ktmp,[2,3]);
        ste_info.md{i}.m=ste_info.md{i}.m(:,:,find(mask_gre_cl));
        ste_info.md{i}.gb0=ste_info.md{i}.gb0(:,find(mask_gre_cl));
        ste_info.md{i}.db0=ste_info.md{i}.db0(find(mask_gre_cl));
    end
    ste_info=setfield(ste_info,'mask_gre',mask_gre);
end