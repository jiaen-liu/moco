function fcrptdm=ste_def_crpt_m(mot_par,par,sorted_ste,res_ste)
    si=size(mot_par);
    if numel(si)>2
        mot_par=reshape(mot_par,[6,prod(si(2:end))]);
        nr=si(2);
        nv=si(3);
    else
        nr=sorted_ste.nparv_cyc;
        nv=si(2)/nr;
    end
    mot_par_fil=gauss_fil_fft(mot_par,nr/2,1,2);
    mot_par_fil=reshape(mot_par_fil,[6,nr,nv]);
    % define corrupted volumes based on motion
    fcrptdm=max(max(mot_par_fil(1:3,:,:),[],2)-min(mot_par_fil(1:3,:,:),[],2),[],1)>=par.max_a | ...
            max((max(mot_par_fil(4:6,:,:),[],2)-min(mot_par_fil(4:6,:,:),[],2)).*res_ste',[],1)>=par.max_t;
    fcrptdm=fcrptdm(:);
    % last volume is not fully sampled
    if length(fcrptdm)>sorted_ste.nvf
        fcrptdm(end)=1;
    end
end