% History
% 2022-12-21: Jiaen Liu: added philips support for reversed PE order for EPI
function [im_recon,flag,epsil,b,k,rr,t]=recon_mb_epi_beta(yn,para,varargin)
% yn is nr by nk (ky*kz) by nch
% some constants
    im_recon=[];
    flag=[];
    epsil=[];
    b=[];
    k=[];
    rr=[];
    st=[];
    t=[];
    nvargin=length(varargin);
    gma=42.58e6*2*pi;
    % nufft interpolation neighbor size
    if numel(para.j)==1
        J=ones(1,para.dimen)*para.j;
    else
        J=para.j;
    end
    % oversampling
    if numel(para.kos)==1
        kos=ones(1,para.dimen,1)*para.kos;
    else
        kos=para.kos(:).';
    end
    % cg iteration max
    n_iter_reset=para.n_iter_reset;
    n_iter=para.n_iter;
    resr=para.resr*1e-3;
    resp=para.resp*1e-3;
    ress=para.ress*1e-3;
    D_res=diag([resr,resp,ress]);
    
    pcaflag = para.pcaflag;
    % 
    nk_shot=para.nk_shot;
    n_interl=para.n_interleaves;
    %
    ie=1;
    irep=1;
    if nvargin > 0
        ie=varargin{1}(1);
        if numel(varargin{1})>1
            irep=varargin{1}(2);
        end
    end
    precision='double';
    if nvargin>1
        precision=varargin{2};
    end
    ro_pol=para.ro_pol;
    ro_pol_ref=para.ro_pol_ref;
    if para.isgre
        if para.n_contrasts>1
            te=para.te_contr(ie)*1e-3;
        else
            te=para.te*1e-3+(ie-1)*para.echo_spacing*1e-6;
        end
    else
        te=para.te_ro-para.te_ro(floor((nk_shot+2)/2))+...
           para.te_contr(ie);
        te=te(:)*1e-3;
        % te: nk_shotxn_interl
        if isfield(para,'pe_order') && strcmp(para.pe_order,'rev_linear')
            % for philips, the pe order can be reversed
            te=te+...
               [n_interl-1:-1:0]*...
               para.echo_spacing*1e-6/n_interl*para.int_te_shift;
        else
            te=te+...
               [0:n_interl-1]*...
               para.echo_spacing*1e-6/n_interl*para.int_te_shift;
        end
    end
    % number of lines per shot
    nline_shot=double(para.np./para.n_interleaves/...
                      para.sense_rate_p);
    % dimension in read-out direction
    nr=double(para.nr);
    % dimension in phase encoding direction
    np=double(para.np);
    % dimension in slice direction
    ns=double(para.n_partitions);
    % total voxels
    n=nr*np*ns;
    % size
    if para.dimen == 2
        si=[nr,np];
    elseif para.dimen==3
        si=[nr,np,ns];
    end
    nch=para.n_channels;
    % total motion status
    nm=para.nm;
    yn=reshape(yn,[numel(yn)/para.n_channels,para.n_channels]);
    % transpose for fast calculation using intel mkl sparse multiplication
    yn=yn.';
    % transpose b1n
    if isfield(para,'b1n') && ~isempty(para.b1n)
        b1n=para.b1n;
        b1n=reshape(b1n,[numel(b1n)/nch/nm,nch,nm]);
        b1n=permute(b1n,[2,1,3]);
        para.b1n=b1n;
        clear b1n;
        b0=para.b0;
        nori_b0=numel(para.b0)/nm;
        b0=reshape(b0,[numel(b0)/nm,nm]);
        if para.isgre
            para.b0=exp(1i*2*pi*b0*te);
        else
            para.b0=exp(1i*2*pi*b0*para.te_contr(ie)*1e-3);
        end
        
        if pcaflag
            npca = para.npca;
        else
            npca=1;
        end
        clear b0;
    end
    % fov unit: meter
    lr=para.fovr*1e-3;
    lp=para.fovp*1e-3;
    ls=para.sthickness*(para.sl_oversamp+1)*1e-3;
    % dwell time: second
    dt=2*para.t_dwell*1e-9;
    %
    pi2i=2*pi*1i;
    pi2=2*pi;
    % 1. create the elements for cg recon based on nufft
    %% prepare nufft parameters based on the motion and field input
    md=para.md;
    nk=0;
    for im=1:nm
        nk=nk+md(im).nk;
    end
    nshot_rep=para.n_interleaves*para.n_partitions/para.sense_rate_s;
    % number of k-lines per repetition
    nk=min(nk,nk_shot*nshot_rep);
    % k space coordinate
    k=zeros(3,nr,nk);
    dp=zeros(nr,nk,precision);
    kr_o=2*pi*((0:nr-1)-floor(nr/2))/lr;
    % kr_o=2*pi*((0:nr-1)-floor(nr/2));
    kr_o_int=(0:nr-1)-floor(nr/2);
    %
    ik=1;
    % 
    for im=1:nm
        % in case of multiple repetitions
        % the number of k-lines is less than one repetition
        for i=1:min(nk,md(im).nk)
            % calculate the new kspace coordinate due to rotation and
            % B0 gradient
            % k space position in clusters is defined
            % in ste_motion.m
            ii=i+nk*(irep-1);
            ishot_m=mod(floor((ii-1)/nk_shot),nshot_rep)+1+...
                    nshot_rep*(irep-1);
            ikline_shot=mod(ii-1,nk_shot)+1;
            if para.isgre
                i_te_interl=1;
                sign_dt=conditional(ro_pol_ref(mod(ie-1,para.nk_shot_ref)+1,...
                                               floor((ie-1)/para.nk_shot_ref)+1),-1,1);
            else
                idx_shot=para.idx_reorder_shot(floor((ik-1)/nk_shot)+1);
                i_te_interl=mod(idx_shot-1,n_interl)+1;
                % Deal with polarity of readout lines
                sign_dt=conditional(ro_pol(ikline_shot,ie),-1,1);
            end
            
            kp_o=md(im).k(1,ii)*pi2/lp;
            ks_o=md(im).k(2,ii)*pi2/ls;
            m=md(im).m(:,:,ishot_m); % 6-parameter rigid motion
            R=rot3d(m(:,1)); % in rad
            dr=m(:,2); % in meter
            ko=[kr_o;ones(1,nr)*kp_o;ones(1,nr)*ks_o];
            gb0=md(im).gb0(:,ishot_m); % gradient of B0 change: unit (Hz/m)
            db0=md(im).db0(ishot_m); % offset of B0 change (Hz)
            dk=(sign_dt*dt*pi2)*(gb0*kr_o_int)+...
               (te(ikline_shot,i_te_interl)*pi2)*gb0;
            % k(:,:,ik)=D_res*R.'*(ko+dk);
            k(:,:,ik)=D_res*(R.'*ko+dk);
            % calculate the phase change due to translation and global
            % B0 change
            dp(:,ik)=exp(1i*((ko+dk).'*dr+...
                             db0*(te(ikline_shot,i_te_interl)+...
                                  sign_dt*dt*kr_o_int.')*pi2));
            ik=ik+1;
        end
    end
    k=-reshape(k,[3,nr*nk]); % -i for nufft
                             % create nufft
    st=nufft_init_efficient(k(1:para.dimen,:).',si,J,...
                            kos.*si,precision,(si-1)/2);
    % st=nufft_init(k(1:para.dimen,:).',si,J,para.kos*si,(si-1)/2);
    % the returned sparse matrix in st.p is in row major mode with 
    % coordinate (coo)
    % convert the sparse matrix coo and its transpose to csr format
    disp('*** Prepare sparse matrix for Intel MKL ***');
    % Estimate memory need for MKL sparse matrix and its transpose
    % column index is int64 8 bytes
    % value array is single or double
    byte_float=4;
    if strcmp(precision,'double')
        byte_float=8;
    end
    mkl_mat_size=0;
    % only count column index+data
    % the rest is minor
    mkl_mat_size=mkl_mat_size+2*st.M*st.nnz_row*(8+byte_float);
    mkl_mat_size=mkl_mat_size/1024^3;
    disp(num2str(memory_linux(),'*** Available memory: %g GB ***'));
    disp(num2str(mkl_mat_size,'*** Needs at least %g Gbyte RAM for MKL matrix ***'));
    % 2024/03/22, Jiaen Liu:
    % mkl_sp_transpose is not memory efficient
    % the program can crash with large image matrix size
    % start segmenting the reconstruction data in case of large image matrix size
    nseg=1;
    if ~isfield(para,'n_seg_mem')
        para.n_seg_mem=1;
    end
    if nm>1
        % if k-means is used, segmentaion follows k-means clustering.
        nseg=1;
    elseif nm==1 && para.n_seg_mem>1
        % sequencial segmentation
        nseg=para.n_seg_mem;
    end
    idx_seg=zeros(2,nseg*nm);
    if nseg==1 
        idx_seg(1,1)=1;
        % deal with multiple repetition
        idx_seg(2,1)=min(para.md(1).nk,nk)*nr;
        if nm>1
            for im=2:nm
                for i=1:im-1
                    idx_seg(1,im)=idx_seg(1,im)+para.md(i).nk*nr;
                end
                idx_seg(2,im)=idx_seg(1,im)+para.md(im).nk*nr;
                idx_seg(1,im)=idx_seg(1,im)+1;
            end
        end
    else
        nk_tot=min(para.md(1).nk,nk)*nr;
        nk_seg=ceil(nk_tot/nr/nseg)*nr;
        for iseg=1:nseg
            idx_seg(1,iseg)=nk_seg*(iseg-1)+1;
            if iseg==nseg
                idx_seg(2,iseg)=nk_tot;
            else
                idx_seg(2,iseg)=nk_seg*iseg;
            end
        end
    end
    st.csr=cell(nseg*nm,1);
    sp_seg=cell(nseg*nm,1);
    nnz_row=st.nnz_row;
    disp(num2str(nseg*nm,'*** Slicing sparse matrix to %d segments! ***'));
    for iseg=1:nseg*nm
        idxtmp1=nnz_row*(idx_seg(1,iseg)-1)+1;
        idxtmp2=nnz_row*idx_seg(2,iseg);
        sp_seg{iseg}.mm=st.p.mm(idxtmp1:idxtmp2);
        sp_seg{iseg}.mm=sp_seg{iseg}.mm-sp_seg{iseg}.mm(1)+1;
        sp_seg{iseg}.kk=st.p.kk(idxtmp1:idxtmp2);
        sp_seg{iseg}.uu=st.p.uu(idxtmp1:idxtmp2);
    end
    % clear p to save memory
    st.p=[];
    ncol=st.N;
    nrow=st.M;    
    % tic;
    for iseg=1:nseg*nm
        % use MKL to transpose sparse matrix
        nrow_seg=sp_seg{iseg}.mm(end);
        st.csr{iseg}.val=sp_seg{iseg}.uu;
        % clear sp_seg{iseg}
        sp_seg{iseg}.uu=[];
        [st.csr{iseg}.col_ind,...
         st.csr{iseg}.row_ptr]=...
            coo2csr(sp_seg{iseg}.kk,...
                    sp_seg{iseg}.mm,...
                    nrow_seg,0,st.nnz_row);
        % clear sp_seg{iseg}
        sp_seg{iseg}.mm=[];
        sp_seg{iseg}.kk=[];
        % only transpose now

        [st.csr{iseg}.val_h,...
         st.csr{iseg}.col_ind_h,...
         st.csr{iseg}.row_ptr_h]=...
            mkl_sp_transpose(conj(st.csr{iseg}.val),...
                             st.csr{iseg}.col_ind,...
                             st.csr{iseg}.row_ptr,nrow_seg,ncol);
        st.csr{iseg}.m=int64(nrow_seg);
        st.csr{iseg}.n=int64(ncol);
        % clear sp_seg{iseg}
        sp_seg{iseg}=[];
    end
    % toc;
    clear sp_seg;
    disp(num2str(memory_linux(),'*** Available memory: %g GB ***'));
    %
    para=setfield(para,'idx_seg',idx_seg);
    para=setfield(para,'st',st);
    clear st;
    para=setfield(para,'dp',dp);
    clear dp;
    %% start iterative reconstruction
    disp 'Starting iterative reconstruction...'
    tic;
    Amult=@(x)recon_mb_epi_cg_beta(x,para,1);

    b=reshape(recon_mb_epi_cg_beta(yn,para,0),[n,1]);
    
    x0=complex(zeros(n,1,precision));

    rv_all=zeros(n_iter*n_iter_reset,1);
    iter_all=0;
    fprintf('The current iteration is (/%d): ',n_iter);
    for i=1:n_iter
        [im_recon,flag,rr,iter,rv]=pcg(Amult,b,1e-3,n_iter_reset,[],[], ...
                                       x0);
        iter_all=iter_all+iter;
        rv_all((i-1)*(n_iter_reset)+1:(i-1)*(n_iter_reset)+iter)=rv(2:iter+1);
        if flag == 0
            break;
        end
        x0=reshape(im_recon,[n,1]);
        fprintf('%d,', i);
        
    end
    t=toc;
    fprintf('\n');
    disp(['*** It took ', num2str(t) ' seconds to finish one reconstruction ***']);
    epsil=rv_all(1:iter_all);
    im_recon=reshape(im_recon,[nr,np,ns]);
end

function y=recon_mb_epi_cg_beta(x,para,retflag)
% solving AT*y=AT*A*x
% retflag: 0, AT*y; 1, AT*A*x
    if isa(x,'single')
        precision='single';
    elseif isa(x,'double')
        precision='double';
    end
    
    pcaflag = para.pcaflag;
    
    if pcaflag
        npca = para.npca;
    else
        npca=1;
    end
    
    b1x=[];
    % dimension in read-out direction
    nr=double(para.nr);
    % dimension in phase encoding direction
    np=double(para.np);
    % dimension in slice direction
    ns=double(para.n_partitions);
    % total voxels
    n=nr*np*ns;
    % oversampling
    if numel(para.kos)==1
        kos=ones(1,para.dimen,1)*para.kos;
    else
        kos=para.kos(:).';
    end
    dimen=para.dimen;
    nos=n*prod(kos);
    % if input x is zero, return zero
    if nnz(x)<1 
        y=zeros(n,1,'like',x);
        return;
    end
    % total motion status
    % or memory segments
    nm=para.nm;
    nseg=1;
    if nm>1
        % if k-means is used, segmentaion follows k-means clustering.
        nseg=1;
    elseif nm==1 && para.n_seg_mem>1
        % sequencial segmentation
        nseg=para.n_seg_mem;
    end
    % fov unit: meter
    lr=para.fovr*1e-3;
    lp=para.fovp*1e-3;
    
    ls=para.sthickness*(para.sl_oversamp+1)*1e-3;
    % dwell time: second
    dt=2*para.t_dwell*1e-9;
    % number of channels
    nch=para.n_channels;
    md=para.md;
    % total lines
    nline=0;
    for im=1:nm
        nline=nline+md(im).nk;
    end
    clear md;
    nline_shot=para.nk_shot;
    nline_rep=para.n_interleaves*para.n_partitions/para.sense_rate_s*...
              nline_shot;
    nline=min(nline,nline_rep);
    %
    pi2i=2*pi*1i;
    pi2=2*pi;
    sn=reshape(para.st.sn,[1,nr,np,ns]);
    idx_seg=para.idx_seg;
    nori_b1=numel(para.b1n)/nch/nm;
    nori_b0=numel(para.b0)/nm;

    if isa(x,'double')
        xm=complex(zeros(1,nr,np,ns));
    else
        xm=complex(zeros(1,nr,np,ns,'single'));
    end
    % Generate index to loop through channels
    % Note: it has to loop through channels because of 
    % memory issue.
    
    % calculate maximum channels to process based on memory usage
    freemem=memory_linux();
    mem_need=0;
    if isa(x,'double')
        mem_need=nr*np*ns*kos(1)*kos(2)*kos(3)*16/1024^3;
    else
        mem_need=nr*np*ns*kos(1)*kos(2)*kos(3)*8/1024^3;
    end
% $$$     size_para=whos('para');
% $$$     size_para=size_para.bytes/1024^3;
% $$$     max_nch_seg_mem=floor((freemem-size_para)*0.6/mem_need);
    max_nch_seg_mem=floor(freemem*0.8/4/mem_need);
    
    if ~isfield(para,'max_nch_seg')
        max_nch_seg=nch;
    else
        max_nch_seg=para.max_nch_seg;
    end
% $$$     disp(['*** Available memory is ',num2str(round(freemem)),' GB ***']);
% $$$     disp(['*** Requested max number of channels is ',num2str(max_nch_seg),' ***']);
    max_nch_seg=min(max_nch_seg,max_nch_seg_mem);
% $$$     disp(['*** Memory supports max number of channels of ',num2str(max_nch_seg),' ***']);
    
    
    n_ch_seg=ceil(nch/max_nch_seg);
    idx_ch_seg=zeros(2,n_ch_seg);
    for i=1:n_ch_seg
        idx_ch_seg(1,i)=(i-1)*max_nch_seg+1;
        if i<n_ch_seg
            idx_ch_seg(2,i)=i*max_nch_seg;
        else
            idx_ch_seg(2,i)=nch;
        end
    end

    if retflag==0
        kd_fwd=conj(para.dp(:).').*reshape(x,[nch,numel(x)/nch]);
    elseif retflag==1
        x=reshape(x,[1,nr,np,ns]);
    end
    
    for im=1:nm
        print_countdown(nm,im,'K-means cluster #:');
        %% prepare the pca data or the k-means B0 and B1 data with nm clusters
        if pcaflag
            %%%%%%%Save scores (PCA Modes) as a list of fields rather than
            %%%%%%%a list of concatonated vectors
            scores = zeros(npca,nr,np,ns);
            for mode = 1:npca
                if isfield(para,'intp_ker_ste') && ~isempty(para.intp_ker_ste)
                    % interpolate scores and coefs
                    score=sparse_csr_mm_prit(para.intp_ker_ste.val,...
                                             para.intp_ker_ste.col_ind,...
                                             para.intp_ker_ste.row_ptr,...
                                             n,...
                                             nori_b0,...
                                             para.scorefilt(:,mode).');
                    score=reshape(score,[nr,np,ns]);

                else
                    score=reshape(para.scorefilt(:,mode),[nr,np,ns]);
                end
                scores(mode,:,:,:) = score;
            end
            coefs = para.coefsinterp;
        else
            if isfield(para,'intp_ker_ste') && ~isempty(para.intp_ker_ste)
                % interpolate b0
                b0=sparse_csr_mm_prit(para.intp_ker_ste.val,...
                                      para.intp_ker_ste.col_ind,...
                                      para.intp_ker_ste.row_ptr,...
                                      n,...
                                      nori_b0,...
                                      para.b0(:,im).');
                b0=reshape(b0,[1,nr,np,ns]);
            else
                b0=reshape(para.b0(:,im),[1,nr,np,ns]);
            end
        end
        
        if isfield(para,'intp_ker_ext') && ~isempty(para.intp_ker_ext)
            % interplate B1
            b1n=sparse_csr_mm_prit(para.intp_ker_ext.val,...
                                   para.intp_ker_ext.col_ind,...
                                   para.intp_ker_ext.row_ptr,...
                                   n,...
                                   nori_b1,...
                                   para.b1n(:,:,im));
            b1n=reshape(b1n,[nch,nr,np,ns]);
        elseif isfield(para,'intp_ker_ste') && ~isempty(para.intp_ker_ste)
            b1n=sparse_csr_mm_prit(para.intp_ker_ste.val,...
                                   para.intp_ker_ste.col_ind,...
                                   para.intp_ker_ste.row_ptr,...
                                   n,...
                                   nori_b1,...
                                   para.b1n(:,:,im));
            b1n=reshape(b1n,[nch,nr,np,ns]);
        else
            b1n=reshape(para.b1n(:,:,im),[nch,nr,np,ns]);
        end
        if retflag==1
            x_sn=x.*sn;
            if pcaflag~=1
                b0_x_sn=b0.*x_sn;
            end
        end
        
        %% first calculate A*x with retflag of 1
        for ich_seg=1:n_ch_seg
            idx_ch=idx_ch_seg(1,ich_seg):idx_ch_seg(2,ich_seg);
            n_ch_seg_cur=numel(idx_ch);
            n_ch_seg_cur=length(idx_ch);
            
            if retflag==1
                if pcaflag
                    kd_fwd=0;
                    % A*x for pca
                    for mode=1:npca
                        kd_fwd_mode=b1n(idx_ch,:,:,:).*(scores(mode,:,:,:).*x_sn);
                        kd_fwd_mode=fft(kd_fwd_mode,[nr]*kos(1),2);
                        kd_fwd_mode=fft(kd_fwd_mode,[np]*kos(2),3);
                        if para.dimen == 3
                            kd_fwd_mode=fft(kd_fwd_mode,[ns]*kos(3),4);
                        end
                        kd_fwd_mode=reshape(kd_fwd_mode,[n_ch_seg_cur,nos]);
                        kd_fwd_tmp=zeros(n_ch_seg_cur,n);
                        for iseg=1:nseg
                            idx_act_seg=iseg+nseg*(im-1);
                            kd_fwd_tmp(:,idx_seg(1,idx_act_seg):idx_seg(2,idx_act_seg))=...
                                sparse_csr_mm_prit(para.st.csr{idx_act_seg}.val,...
                                                   para.st.csr{idx_act_seg}.col_ind,...
                                                   para.st.csr{idx_act_seg}.row_ptr,...
                                                   para.st.csr{idx_act_seg}.m,...
                                                   para.st.csr{idx_act_seg}.n,...
                                                   kd_fwd_mode);
                        end
                        kd_fwd_mode=kd_fwd_tmp;
                        clear kd_fwd_tmp;
                        kd_fwd_mode=reshape(coefs(:,mode),[1,1,nline]).*...
                            reshape(kd_fwd_mode,[n_ch_seg_cur,nr,nline]);
                        kd_fwd=kd_fwd+kd_fwd_mode;
                    end
                else
                    % A*x for non-pca 
                    kd_fwd=b1n(idx_ch,:,:,:).*(b0_x_sn);
                    kd_fwd=fft(kd_fwd,[nr]*kos(1),2);
                    kd_fwd=fft(kd_fwd,[np]*kos(2),3);
                    if para.dimen == 3
                        kd_fwd=fft(kd_fwd,[ns]*kos(3),4);
                    end
                    kd_fwd=reshape(kd_fwd,[n_ch_seg_cur,nos]);
                end
            end % if retflag == 1
            % loop over k-space segments and pca components
            % for the second half of the forward calculation
            for iseg=1:nseg
                idx_act_seg=iseg+nseg*(im-1);
                nkline_seg=(idx_seg(2,idx_act_seg)-idx_seg(1,idx_act_seg)+1)/nr;
                idx_seg_line_start=(idx_seg(1,idx_act_seg)-1)/nr+1;
                idx_seg_line_end=idx_seg(2,idx_act_seg)/nr;
                for mode=1:npca
                    if pcaflag==1
                        kd=reshape(conj(coefs(idx_seg_line_start:idx_seg_line_end,mode)),...
                                   [1,1,nkline_seg]).*...
                           reshape(kd_fwd(idx_ch,idx_seg(1,idx_act_seg):idx_seg(2,idx_act_seg)),...
                                   [n_ch_seg_cur,nr,nkline_seg]);
                        kd=reshape(kd,[n_ch_seg_cur,nr*nkline_seg]);
                    else
                        if retflag==1
                            kd=sparse_csr_mm_prit(para.st.csr{idx_act_seg}.val,...
                                                  para.st.csr{idx_act_seg}.col_ind,...
                                                  para.st.csr{idx_act_seg}.row_ptr,...
                                                  para.st.csr{idx_act_seg}.m,...
                                                  para.st.csr{idx_act_seg}.n,...
                                                  kd_fwd);
                        else
                            kd=kd_fwd(idx_ch,idx_seg(1,idx_act_seg):idx_seg(2,idx_act_seg));
                        end
                    end
                    
                    
                    
                    kd=sparse_csr_mm_prit(para.st.csr{idx_act_seg}.val_h,...
                                          para.st.csr{idx_act_seg}.col_ind_h,...
                                          para.st.csr{idx_act_seg}.row_ptr_h,...
                                          para.st.csr{idx_act_seg}.n,...
                                          para.st.csr{idx_act_seg}.m,...
                                          kd);
                    kd=reshape(kd,[n_ch_seg_cur,...
                                   nr*kos(1),...
                                   np*kos(2),...
                                   ns*kos(3)]);
                    kd=ifft(kd,nr*kos(1),2);
                    kd=ifft(kd,np*kos(2),3);
                    if dimen == 3
                        kd=ifft(kd,ns*kos(3),4);
                    end
                    kd=kd(:,1:nr,1:np,:);
                    if dimen==3
                        kd=kd(:,:,:,1:ns);
                    end
                    if pcaflag
                        xm=xm+...
                           sum(kd.*conj(b1n(idx_ch,:,:,:)),1).*...
                           conj(scores(mode,:,:,:));
                    else
                        xm=xm+...
                           sum(kd.*conj(b1n(idx_ch,:,:,:)),1).*...
                           conj(b0);
                    end
                end % loop over pca
            end % loop over memory segments
        end % loop over ch seg
    end % loop over k-means
    xm=(xm.*conj(sn))*nos;
    y=reshape(xm,[n,1]);
end
