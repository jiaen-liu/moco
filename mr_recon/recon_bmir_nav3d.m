function [dnav,pe,im_nav]=recon_bmir_nav3d(mid,varargin)
    p=inputParser;
    no_fov_crct=0;
    no_pc=0;
    addParameter(p,'no_fov_crct',no_fov_crct,@isnumeric);
    addParameter(p,'no_pc',no_pc,@isnumeric);
    p.parse(varargin{:});
    no_fov_crct=p.Results.no_fov_crct;
    no_pc=p.Results.no_pc;
    para=extract_para_philips(mid);
    dnav=[];
    pe=[];
    if ~para.b_ste_en
        return;
    end
    dy=[0,2,1,3];
    dz=[0,1];
    s1nav=para.n_interleaves_steref;
    s2nav=2;
    r=s1nav*s2nav;
    nch=para.n_channels;
    nrnav=para.steref_dim_r*2;
    npnav=para.steref_dim_p;
    nsnav=para.steref_dim_s;
    nshot=para.n_main_tr;
    nechonav=npnav/para.n_interleaves_steref;
    n_shot_fast_nav3d=nsnav/s2nav+1;
    n_shot_full_nav3d=n_shot_fast_nav3d*r;
    % define phase encoding steps
    py=zeros(n_shot_fast_nav3d*nechonav,r);
    pz=zeros(n_shot_fast_nav3d*nechonav,r);
    for i=1:r
        pe=gen_pe_sense(npnav,nsnav,s1nav,s2nav,1,dy(mod(i-1,s1nav)+1),dz(floor((i-1)/s1nav)+1));
        py(nechonav+1:end,i)=pe(1,:);
        pz(nechonav+1:end,i)=pe(2,:);
    end
    py=[py(:);zeros(nshot*nechonav-n_shot_full_nav3d*nechonav,1)];
    pz=[pz(:);zeros(nshot*nechonav-n_shot_full_nav3d*nechonav,1)];
    for i=1+n_shot_full_nav3d*nechonav:nshot*nechonav
        py(i)=py(mod(i-1,n_shot_fast_nav3d*nechonav)+1);
        pz(i)=pz(mod(i-1,n_shot_fast_nav3d*nechonav)+1);
    end
    py=py(:);
    pz=pz(:);
    pe=[py.';pz.'];
    [dnav,hnav]=read_raw_philips(mid,'mix',1,'type',1);
    dnav=reshape(dnav,[nrnav,nch,nechonav,nshot]);
    % reformat to nrnav,nechonav,nch,nshot
    dnav=permute(dnav,[1,3,2,4]);
    % truncate
    dnav=fftmr(dnav,-1,1);
    dnav=dnav(idx_truncate(nrnav,nrnav/2),:,:,:);
    nrnav=nrnav/2;
    % phase correction
    if ~no_pc
        if abs(para.frequency/42.58e6-7)<0.5
            %% 7T data perform phase correction in two receive channel groups
            idx=[1:8:nch]+[0:3].';
            idx1=idx(:);
            idx=[5:8:nch]+[0:3].';
            idx2=idx(:);
            idx=[idx1,idx2];
            necho_pc=nechonav;
            pc_by_echo=1;
        else
            idx=[1:nch].';
            necho_pc=1;
            pc_by_echo=0;
        end
        for ishot=1:nshot
            if mod(ishot-1,n_shot_fast_nav3d)==0
                pha=zeros(nrnav,necho_pc,size(idx,2));
                if pc_by_echo
                    for i=1:size(idx,2)
                        for iepc=1:necho_pc
                            if iepc==1
                                dblpls = squeeze(dnav(:,1:3,idx(:,i),ishot));
                                pha(:,iepc,i)=dpOddEven(dblpls,2);
                            elseif iepc==nechonav
                                dblpls = squeeze(dnav(:,nechonav-2:nechonav,idx(:,i),ishot));
                                pha(:,iepc,i)=dpOddEven(dblpls,2);
                                if mod(nechonav-2,2)==0
                                    pha(:,iepc,i)=-1*pha(:,iepc,i);
                                end
                            else
                                dblpls = squeeze(dnav(:,iepc-1:iepc+1,idx(:,i),ishot));
                                pha(:,iepc,i)=dpOddEven(dblpls,2);
                                if mod(iepc-1,2)==0
                                    pha(:,iepc,i)=-1*pha(:,iepc,i);
                                end
                            end
                        end
                    end
                else
                    for i=1:size(idx,2)
                        dblpls = squeeze(dnav(:,:,idx(:,i),ishot));
                        pha(:,1,i)=dpOddEven(dblpls,2);
                    end
                end
            else
                for i=1:size(idx,2)
                    if pc_by_echo
                        dnav(:,1:2:end,idx(:,i),ishot)=dnav(:,1:2:end,idx(:,i),ishot)./...
                            exp(1i*pha(:,1:2:end,i)/2);
                        dnav(:,2:2:end,idx(:,i),ishot)=dnav(:,2:2:end,idx(:,i),ishot).*...
                            exp(1i*pha(:,2:2:end,i)/2);                        
                    else
                        dnav(:,1:2:end,idx(:,i),ishot)=dnav(:,1:2:end,idx(:,i),ishot)./...
                            exp(1i*pha(:,i)/2);
                        dnav(:,2:2:end,idx(:,i),ishot)=dnav(:,2:2:end,idx(:,i),ishot).*...
                            exp(1i*pha(:,i)/2);
                    end
                    
                end
            end
        end
    end
    %% test code
% $$$     if abs(para.frequency/42.58e6-7)<0.5
% $$$         dnav(:,2:2:end,:,:)=-dnav(:,2:2:end,:,:);
% $$$     end
    dnav=fftmr(dnav,1,1);
    if ~no_fov_crct
        % correct fov
        % readout
        lin_pha_ro=([0:nrnav-1].'-floor(nrnav/2))*2*pi*para.m_shift(1)/para.fovr;
        dnav=dnav.*exp(-1i*lin_pha_ro);
        % phase direction
        dnav=dnav.*exp(-1i*2*pi*reshape(py,[1,nechonav,1,nshot])*para.p_shift(1)/para.fovp);
        % slice direction
        dnav=dnav.*exp(-1i*2*pi*reshape(pz,[1,nechonav,1,nshot])*para.s_shift(1)/para.fovs);
    end
    % create a image from the first fully sampled volume
    im_nav=zeros(nrnav,npnav,nsnav,nch);
    for i=1:n_shot_full_nav3d
        if mod(i-1,n_shot_fast_nav3d)==0
            continue;
        end
        for j=1:nechonav
            idx=j+(i-1)*nechonav;
% $$$         if py(idx)==0 && pz(idx)==0
% $$$             continue;
% $$$         end
            m(mod(py(idx)+npnav/2,npnav)+1,mod(pz(idx)+nsnav/2,nsnav)+1)=1;
            im_nav(:,mod(py(idx)+npnav/2,npnav)+1,...
                   mod(pz(idx)+nsnav/2,nsnav)+1,...
                   :)=dnav(:,j,:,i);
        end
    end
    im_nav=fftmr(im_nav,-1,[1,2,3]);
end