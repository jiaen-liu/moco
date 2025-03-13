% 2022-12-23: Jiaen Liu, added support for Philips
function im=proc_sorted_ste(sorted_ste,varargin)
    fullfov=0;
    no_comb=0;
    p=inputParser;
    b1=[];
    icontr=[];
    gcalib=[];
    idx_f=[];
    apodization=0.25;
    vendor='siemens';
    addParameter(p,'fullfov',fullfov,@isnumeric);
    addParameter(p,'no_comb',no_comb,@isnumeric);
    addParameter(p,'sense',b1,@isnumeric);
    addParameter(p,'icontr',icontr,@isnumeric);
    addParameter(p,'grappa',gcalib,@isstruct);
    addParameter(p,'idx_f',idx_f,@isnumeric);
    addParameter(p,'apodization',apodization,@isnumeric);
    addParameter(p,'vendor',vendor,@ischar);
    p.parse(varargin{:});
   
    fullfov=p.Results.fullfov;
    no_comb=p.Results.no_comb;
    b1=double(p.Results.sense);
    icontr=p.Results.icontr;
    gcalib=p.Results.grappa;
    idx_f=p.Results.idx_f;
    apodization=p.Results.apodization;
    vendor=p.Results.vendor;
    
    if isfield(sorted_ste,'vendor')
        % bypass input is vendor is defined in preprocessed data
        vendor=sorted_ste.vendor;
        if iscell(vendor)
            vendor=vendor{1};
        end
    end
    
    if ~isempty(gcalib)
        w=gcalib.w;
        parg=gcalib.par;
    end
    % parellel imaging parameters
    getvar_struct(sorted_ste,'s1','s2','necho',...
                             'nshot','inv_cov',...
                             'ncontr','dkz','d_nblpl','k_nblpl',...
                             'nshot_nblpl','nv','nvf',...
                             'nparv_cyc');
    
    nr=sorted_ste.dims_k(1);
    np=sorted_ste.dims_k(2);
    ns=sorted_ste.dims_k(3);
    nch=sorted_ste.dims_k(4);
    
    n=nr*np*ns;
    n_nch=n*nch;

    % number of echoes per contrast
    necho_contr=necho/ncontr;
    % r is the number of partial volumes to cycle through a full
    % fov; r may depend on the ste parameters.
    r=nparv_cyc;
    % fullfov is not 1 or 0 then fullfov is used to decrease the accelleration ratio
    if fullfov == 2 && sorted_ste.y_cyc ~= 2
        error('The acquisision needs to be interleaved to combine data');
    end
    rpe=r;
    s1pe=s1;
    if fullfov ~= 1 && fullfov ~= 0
        r=r/fullfov;
        s1=s1/fullfov;
    end
    %
    if fullfov~=0
        nshot_v=ns/s2*fullfov;
    else
        nshot_v=ns/s2;
    end
    nshot_fv=nshot_v*r;
    % caipi reconstruct matrix
    smat=cell(r,1);

    if fullfov ~= 1 && isempty(b1) && isempty(gcalib)
        error('Provide a sensitivity map');
    end
    % SENSE inversion matrix
    if fullfov ~= 1
        if sorted_ste.y_cyc==1
            dy=mod([0:rpe-1],s1pe);
        elseif sorted_ste.y_cyc == 2
            hsp = floor(s1pe/2);
            dy=[(0:hsp-1);(hsp:2*hsp-1)];
            dy=dy(:);
            if mod(s1pe,2) == 1
                dy=[dy;s1pe-1];
            end
        elseif sorted_ste.y_cyc == 0
            dy=zeros(rpe,1);
        end
        if fullfov~=0 && sorted_ste.y_cyc==1
            dy=dy(1:fullfov);
        elseif fullfov~=0 && sorted_ste.y_cyc==2
            dy=dy(1:fullfov:end);
        end
        if isempty(gcalib)
            if strcmp(vendor,'philips')
                error('*** SENSE recon for the navigator is not supported for Philips data ***');
            end
            for i=1:r
                dz=floor((i-1)/s1);
                smat{i}=sense_mat(b1,inv_cov,s1,s2,dkz,dy(mod(i-1,s1)+1),dz);
            end
        end
    end
    % reconstruct the images
    if isempty(icontr)
        icontr=1:ncontr;
        ncontr_recon=ncontr;
    else
        ncontr_recon=numel(icontr);
    end
    % reconstruct selected segments
    if isempty(idx_f)
        nvf_recon=ceil(nshot_nblpl/nshot_v/r);
        nv=floor(nshot_nblpl/nshot_v);
        iv_recon=[1:nv];
    else
        nvf_recon=length(idx_f);
        nv=nvf_recon*r;
        iv_recon=[1:r].'+(idx_f(:).'-1)*r;
        iv_recon=iv_recon(:);
    end
    if fullfov == 1
        im=zeros(nr,np,ns,nch,ncontr_recon,nvf_recon);
        di=zeros(nr,np,ns,nch,ncontr_recon,nvf_recon);
    else

        if isempty(gcalib)
            % sense
            im=zeros(nr,np,ns,ncontr_recon,r,nvf_recon);
            di=zeros(n,ncontr_recon,r,nvf_recon);
        else
            % grappa
            im=zeros(nr,np,ns,nch,ncontr_recon,r,nvf_recon,'single');
            di=zeros(nr,np*ns/r,nch,ncontr_recon,r,nvf_recon,'single');
            kyz=reshape(k_nblpl,[2,necho_contr,ncontr_recon,nshot_v,r]);
            kyz=kyz(:,:,icontr,:,:);
            kyz=permute(kyz,[1,2,4,3,5]);
            kyz=reshape(kyz,[2,necho_contr*nshot_v,ncontr_recon,r]);
        end
    end
    % assemble data 
    % for full fov, fft is also performed
    % but not for sense or grappa
    for ic=icontr
        for iv=1:nv
            ir=mod(iv-1,r)+1;
            iv_full=floor((iv-1)/r)+1;
            if strcmp(vendor,'siemens')
                kcur=k_nblpl(:,(ic-1)*necho_contr+1:ic*necho_contr, ...
                             (ir-1)*nshot_v+1:ir*nshot_v);
            elseif strcmp(vendor,'philips')
                % The first full volume is the same as the siemens version
                % later on, only accelerated navigator is acquired
                if iv<=r
                    kcur=k_nblpl(:,(ic-1)*necho_contr+1:ic*necho_contr, ...
                                 (ir-1)*nshot_v+1:ir*nshot_v);
                else
                    ir_tmp=1;
                    kcur=k_nblpl(:,(ic-1)*necho_contr+1:ic*necho_contr, ...
                                 (ir_tmp-1)*nshot_v+1:ir_tmp*nshot_v);
                end
            end
            dcur=d_nblpl(:,(ic-1)*necho_contr+1:ic*necho_contr,...
                         :,(iv_recon(iv)-1)*nshot_v+1:iv_recon(iv)*nshot_v);
            % apodize based on the k coordinate
            if apodization~=-1
                k1=mod(floor(np/2)+squeeze(kcur(1,:,:)),np);
                k2=mod(floor(ns/2)+squeeze(kcur(2,:,:)),ns);
                if ns > 1
                    wapd=win_tukey(numel(k1),apodization,k1/(np-1)).*...
                         win_tukey(numel(k2),apodization,k2/(ns-1));
                else
                    wapd=win_tukey(numel(k1),apodization,k1/(np-1));
                end
                wapd=reshape(wapd,[1,necho_contr,1,nshot_v]);
                dcur=dcur.*wapd;
            end
            
            if isempty(gcalib)
                % use sense or full-fov
                if fullfov==1 && ir==1
                    dk=zeros(nr,np,ns,nch);
                elseif fullfov==0
                    dk=zeros(nr,np,ns,nch);
                end
                % store data
                for ishot=1:nshot_v
                    for ie=1:necho_contr
                        dk(:,floor(np/2)+1+kcur(1,ie,ishot),...
                           mod(floor(ns/2)+kcur(2,ie,ishot),ns)+1,:)=...
                            dcur(:,ie,:,ishot);
                    end
                end
                
            else
                % use grappa
                dcur=permute(dcur,[1,2,4,3]);
                di(:,:,:,ic,ir,iv_full)=reshape(dcur,[nr,necho_contr*nshot_v,nch]);
            end

            if fullfov ~= 1
                if isempty(gcalib)
                    % use sense
                    di_tmp=fftmr(dk,-1,[1,2,3])*(nr*np*ns);
                    di_tmp=reshape(di_tmp,[n,nch]);
                    di_tmp=di_tmp*conj(inv_cov);
                    di_tmp=reshape(di_tmp,[nr,np,ns,nch]);
                    di(:,ic,ir,iv_full)=sum(reshape(conj(b1).*di_tmp, ...
                                                        [n,nch]),2);
                end
            elseif ir==r
                di_tmp=fftmr(dk,-1,[1,2,3])*(nr*np*ns);
                im(:,:,:,:,ic,iv_full)=di_tmp;
            end
            % im(:,:,:,ic,iv)=reshape(smat{ir}\di,[nr,np,ns]);
        end
    end
    % reconstruction
    if fullfov ~= 1
        if isempty(gcalib)
            % sense
            for ir=1:r
                disp(['Partitial Volume #:' int2str(ir)]);
                im(:,:,:,:,ir,:)=reshape(smat{ir}\reshape(di(:,:,ir,:),[n,ncontr_recon*nvf_recon]),[nr,np,ns,ncontr_recon,nvf_recon]);
            end
        else
            % grappa
            % should segment data to avoid memory overflow
            nv_seg=40;
            nseg=ceil(nvf_recon/nv_seg);
            idx_seg=zeros(2,nseg);
            for iseg=1:nseg
                idx_seg(1,iseg)=nv_seg*(iseg-1)+1;
                idx_seg(2,iseg)=min(nv_seg*iseg,nvf_recon);
            end
            di_cell=cell(nseg);
            for iseg=1:nseg
                di_cell{iseg}=...
                    single(di(:,:,:,:,:,...
                              idx_seg(1,iseg):idx_seg(2,iseg)));
            end

            clear di;
            for iseg=1:nseg
                print_countdown(nseg,iseg,'GRAPPA Segments: ');
                ditmp=di_cell{iseg};
                nv_seg_cur=idx_seg(2,iseg)-idx_seg(1,iseg)+1;
                imtmp=zeros(nr,np,ns,nch,ncontr_recon,r,nv_seg_cur,'single');
                parfor ir=1:r
                    ditmp_parfor=ditmp(:,:,:,:,ir,:);
                    imtmp_parfor=zeros(nr,np,ns,nch,ncontr_recon,nv_seg_cur,'single');
                    kyztmp=[];
                    if strcmp(vendor,'siemens')
                        kyztmp=kyz(:,:,:,ir);
                    elseif strcmp(vendor,'philips')
                        kyztmp=kyz(:,:,:,1);
                    end
                    for iv=1:nv_seg_cur
                        for ic=icontr
                            imtmp_parfor(:,:,:,:,ic,iv)=...
                                grappa_recon(ditmp_parfor(:,:,:,ic,iv),...
                                                              kyztmp(:,:,ic),w,parg);
                        end
                    end
                    imtmp(:,:,:,:,:,ir,:)=imtmp_parfor;
                end
                im(:,:,:,:,:,:,idx_seg(1,iseg):idx_seg(2,iseg))=...
                    imtmp;
            end
            if strcmp(vendor,'philips')
                % reconstructe the first volume with the correct k-space info
                ditmp=di_cell{1}(:,:,:,:,:,1);
                imtmp=zeros(nr,np,ns,nch,ncontr_recon,r,'single');
                for ir=1:r
                    kyztmp=kyz(:,:,:,ir);
                    for ic=icontr
                        imtmp(:,:,:,:,ic,ir)=...
                            grappa_recon(ditmp(:,:,:,ic,ir),...
                                         kyztmp(:,:,ic),w,parg);
                    end
                end
                im(:,:,:,:,:,:,1)=...
                    imtmp;
            end
            clear di_cell imtmp ditmp;
        end
    elseif ~no_comb
        im_mag=sum(abs(im).^2,4).^0.5;
        % im_pha=squeeze(angle(sum(im.*conj(im(:,:,:,:,1,:)),4)));
        im_pha=angle(sum(im.*abs(im),4));
        im=im_mag.*exp(1i*im_pha);
        im=reshape(im,[nr,np,ns,ncontr_recon,nvf_recon]);
    end
end
