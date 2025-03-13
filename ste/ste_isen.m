function isen=ste_isen(freq_fluct,sorted_ste,par)
    nmin=10;
    nvf=sorted_ste.nvf;
    nvf_ceil=length(freq_fluct);
    fcrptd=freq_fluct>par.max_f;
    if nvf_ceil>nvf
        fcrptd(nvf_ceil)=1;
    end
    idx_uncrpt=[1:nvf];
    % uncorrupted
    idx_uncrpt=idx_uncrpt(~fcrptd(1:nvf));
    nuncrpt=length(idx_uncrpt);
    if nuncrpt==0
        error('*** All full fovs are corrupted in the strict criteria! ***');
    end
    % determine the full fovs in these three which covers the k-space
    % closest to the center
    % choose that one
    % determine shots with central kyz
    % kmin=col(min(sum(sorted_ste.kyz(:,1:sorted_ste.para.nk_shot,:).^2,1).^0.5,[],2));
% $$$     [~,idx_kyz_0]=...
% $$$         min(col(min(sum(sorted_ste.kyz(:,1:sorted_ste.para.nk_shot,:).^2,1).^0.5,[],2)));
    % get all k coordinates for shots
    % this may change for non-standard cartesian trajectories
    %%----------------%%

    if isfield(sorted_ste,'kyz')
        ky=sorted_ste.kyz(1,1:sorted_ste.para.nk_shot,:);
        kz=sorted_ste.kyz(2,1:sorted_ste.para.nk_shot,:);
    elseif sorted_ste.para.isgre
        s2=sorted_ste.para.sense_rate_s;
        s1=sorted_ste.para.sense_rate_p;
        np=sorted_ste.para.np;
        n_partitions=sorted_ste.para.n_partitions;
        ky=[0:s1:np-1].'-floor(np/2);
        kz=[0:s2:n_partitions-1].'-floor(n_partitions/2);
        [ky,kz]=ndgrid(ky,kz);
        if sorted_ste.para.dkz_caipi~=0
            error('*** CAIPI is not supported in old data format! ***');
        end
    end
    %%----------------%%
    km=zeros(nuncrpt,1);
    for i=1:nuncrpt
        tmp=find(sorted_ste.idx_shot(1,:)==idx_uncrpt(i));
        km(i)=min(col((ky(1,:,tmp).^2+kz(1,:,tmp).^2).^0.5));
        % km(i)=min((ky(tmp).^2+kz(tmp).^2).^0.5);
    end
    [~,idx_kmin]=mink(km,nmin);
    idx_candi=idx_uncrpt(idx_kmin);
    [~,I]=min(freq_fluct(idx_candi));
    isen=idx_candi(I);
end