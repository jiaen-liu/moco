function mot_par=ste_rect_mot_fast(mot_par,par)
% $$$     mot_par_uc=mot_par(:,:,~fcrptd);
% $$$     mot_par_uc_ave=mean(mot_par_uc,2);
% $$$     mot_par_uc_dif=mean(mot_par_uc-mot_par_uc_ave,3);
% $$$     mot_par=mot_par-mot_par_uc_dif;
    if par.rect_cyc==0
        mot_par_uc=mot_par;
        mot_par_uc_ave=mean(mot_par_uc,2);
        mot_par_uc_dif=median(mot_par_uc-mot_par_uc_ave,3);
        mot_par=mot_par-mot_par_uc_dif;
    elseif par.rect_cyc==1
        si=size(mot_par);
        cyc=si(2);
        mot_par=combine_dim(mot_par,[2,3]);
        mask=false(si(2)*si(3),1);
        mask(si(2)+1:si(2)*(si(3)-1))=true;
        mot_par=steRegress(mot_par,cyc,cyc/2,mask);
        mot_par=reshape(mot_par,si);
    end
end