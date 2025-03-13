function mot_par_comb=ste_comb_motion_beta(ste_info,mot_par_f,mot_par,par,sorted_ste)
% estimate time resolved motion estimation combining artifact-free
% full-fov images with high-temporal resolution accelerated images
    fcrptd=ste_info.fcrptd;
    nv_cyc=sorted_ste.nparv_cyc;
    mot_par_comb=zeros(size(mot_par));
    nv=size(mot_par_f,2);
    if ~isfield(par,'mot_cmb') || par.mot_cmb==0
        % use fast navigator
        mot_par_comb=mot_par;
    elseif par.mot_cmb==2
        % use full navigator
        mot_par_comb=repmat(reshape(mot_par_f,[6,1,nv]),[1,nv_cyc,1]);
    elseif par.mot_cmb==1
        % combine
        mot_par_comb(:,:,~fcrptd)=repmat(reshape(mot_par_f(:,~fcrptd),[6,1,total(~fcrptd)]),[1,nv_cyc,1]);
        mot_par_comb(:,:,fcrptd)=mot_par(:,:,fcrptd);
    end
    % the first volume is not stable
    mot_par_comb(:,:,1)=mot_par_comb(:,:,2);
end