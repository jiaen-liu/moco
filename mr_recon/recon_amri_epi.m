% Modifiction history:
%    2022-09-30..10-03 PvG & JAdZ
%        Fix for data in which mdh only contains positive k-space line numbers.
%    2022-10-04 JAdZ
%        Now calling sort_siemens with the -filemode 0666 option to allow all to
%        read and overwrite.

function [y,para]=recon_amri_epi(mid,varargin)
% check if the data is available
    fd=get_file_filter('.',['MID*',num2str(mid),'.raw.svd']);
    if isempty(fd)
        cmd=['sort_siemens -filemode 0666 ' num2str(mid)];
        if system(cmd)~=0
            error(['*** ',cmd,' was not successful! ***']);
        end
        pause(5);
    end
    p=inputParser;
    k_return=0;
    no_b0_crct=0;
    no_comb=0;
    no_pc=0;
    no_nav_regress=1;
    no_fov_crct=0;
    ord_pha_crct=2;
    mid_b1=[];
    mid_blpo=mid;
    ste=[];
    addParameter(p,'k_return',k_return,@isnumeric);
    addParameter(p,'no_b0_crct',no_b0_crct,@isnumeric);
    addParameter(p,'no_comb',no_comb,@isnumeric);
    addParameter(p,'no_pc',no_pc,@isnumeric);
    addParameter(p,'no_nav_regress',no_nav_regress,@isnumeric);
    addParameter(p,'no_fov_crct',no_fov_crct,@isnumeric);
    addParameter(p,'ord_pha_crct',ord_pha_crct,@isnumeric);
    addParameter(p,'mid_b1',mid_b1,@isnumeric);
    addParameter(p,'mid_blpo',mid_blpo,@isnumeric);
    addParameter(p,'ste',ste,@isstruct);
    p.parse(varargin{:});
    k_return=p.Results.k_return;
    no_b0_crct=p.Results.no_b0_crct;
    no_comb=p.Results.no_comb;
    no_pc=p.Results.no_pc;
    no_nav_regress=p.Results.no_nav_regress;
    no_fov_crct=p.Results.no_fov_crct;
    ord_pha_crct=p.Results.ord_pha_crct;
    mid_b1=p.Results.mid_b1;
    mid_blpo=p.Results.mid_blpo;
    ste=p.Results.ste;
    if k_return
        no_comb=1;
    end
    % define parameters
    acqType='main';
    apodiz=0.25;
    para=extract_para(mid);
    nr_os=para.nr_os;
    nr=para.nr;
    np=para.np;
    ns=para.n_slices;
    ninterl=para.n_interleaves;
    s1=para.sense_rate_p;
    s2=para.sense_rate_s;
    nch=para.n_channels;
    npar=para.n_partitions;
    n_delays=para.n_delays;
    if para.ir_bunch_slices && para.bunch_no_cyc
        n_delays=1;
    end
    nreps=para.n_reps;
    n_echo_contr=para.np/para.n_interleaves/...
        para.sense_rate_p;
    n_echo_shot=0;
    ncontr=para.n_contrasts;
    for i=1:ncontr
        n_echo_shot=n_echo_shot+n_echo_contr+...
            para.n_refs(i);
    end
    if para.isgre
        n_image_echo_shot=n_echo_shot;
    else
        n_image_echo_shot=n_echo_shot-para.n_refs(1);
    end
    n_echo_intl=para.np/para.n_interleaves/s1;
    if contains(para.idea_v,'VB','IgnoreCase',true)
        k=genk(para);
        sik=size(k);
        k=reshape(k,[sik(1),sik(2),1,sik(3)]);
        k=repmat(k,[1,1,para.n_slices,1]);
    else
        k=get_pe_mdh(mid,para);
        sik=size(k);
        k=reshape(k,[sik(1),sik(2)*sik(3),numel(k)/sik(1)/sik(2)/sik(3)]);
    end
	% make sure that this is not a positive-only numbered version (e.g. AMRI_epi ~v1.80)
	if min(k(1,:,:))>=0
		k(1,:,:)=k(1,:,:)-floor(para.np/2);
	end
	if min(k(2,:,:))>=0
		k(2,:,:)=k(2,:,:)-floor(para.n_partitions/2);
	end
    n_echo_nav=0;
    if para.b_nav_en
        n_echo_nav = conditional(para.nav_type==0,1,para.nav_type);
    end
    n_echo_steref=para.n_echo_steref;
    % define file names
    fn_ste=rp(['mid',num2str(mid),'.steref4recon.svd']);
    if ~isempty(mid_b1)
        fn_b1=fullfile(pwd,['mid',num2str(mid_b1),'.b1info.mat']);
    end
    % read out polarity
    [dmdh,mdhTempl]=readmdh(mid,'noise');
    ro_pol=zeros(n_echo_shot,1);
    te=zeros(n_echo_shot,1);
    for i=1:n_echo_shot
        mdhtmp=cast2struct(dmdh(:,1,n_echo_nav+...
                                n_echo_steref+i),...
                           mdhTempl);
        ro_pol(i)=evalmaskbit(mdhtmp,25);
        te(i)=mdhtmp.te;
    end
    % calculate covariance matrix
    cov_mat=covSiem(mid);
    inv_cov=inv(cov_mat);
    % get image data
    % dim: [nx,necho_contr*ncontr,nslice,nch,ntr]
    dk=readSortSiem(mid,acqType,1);
    if n_echo_shot~=size(dk,2)
        error(['*** The number of echos per shot',...
               ' is not consistant! ***']);
    end
    % correct polarity
    dk(:,find(ro_pol==0),:,:,:)=...
        flipdim(dk(:,...
                   find(ro_pol==0),...
                   :,:,:),1);
    % regridding and de-oversampling in readout
    if para.ramp_samp_frac > 0.01
        disp('*** Reconstructing ramp-sampled data ***');
        dx=regridding_arr(dk,para,apodiz);
        clear dk;
    else
        dk=apodize_arr(dk,apodiz,1);
        dx=fftmr(dk,-1,1);
    end
    % nr x necho x nslice x nch x np
    dx=dx(idx_truncate(nr_os,nr),:,:,:,:);
    % correct epi ghost
    % retrieve blipoff data for epi ghost correction
    dblpo = [];
    if para.n_blipoff_reps>0 || mid_blpo~=mid
        disp('*** Retrieving Blipoff data ***');
        dblpo=get_blipoff_data(mid_blpo);
        if mid_blpo~=mid
            disp('*** The blipoff is from another scan ***');
            dblpo=dblpo(:,1:n_echo_contr,:,:);
        end
        para_blpo=extract_para(mid_blpo);
        ro_pol_blpo=para_blpo.ro_pol(:);
% $$$         dblpo=readSortSiem(mid,'blipoff',1);
% $$$         % correct polarity
% $$$         dblpo(:,find(ro_pol==0),:,:,:)=...
% $$$             flipdim(dblpo(:,...
% $$$                           find(ro_pol==0),...
% $$$                           :,:,:),1);
% $$$         % regridding
% $$$         dblpo=regridding_arr(dblpo,para,apodiz);
% $$$         dblpo=dblpo(idx_truncate(nr_os,nr),:,:,:,:);
% $$$         dblpo=reshape(dblpo,[nr,n_echo_shot,ns,...
% $$$                             nch,npar/s2,para.n_blipoff_reps]);
% $$$         % conver to image domain
% $$$         dblpo=fftmr(dblpo,-1,5);
% $$$         dblpo=permute(dblpo,[1,2,4,5,3,6]);
% $$$         % Based on the obersevation that 
% $$$         % eddy current is mostly along read-out direction
% $$$         % other dimentions will be averaged in dpOddEven.m
% $$$         dblpo=mean(dblpo,6);
% $$$         dblpo=reshape(dblpo,[nr,n_echo_shot,...
% $$$                             nch,ns*npar/s2]);
    end
    
    % strip reference scans for EPI ghost correction
    if ~para.isgre && para.n_refs(1)>0
        if para.n_contrasts>1
            error('*** This code does not process multi-contrast data with reference line');
        end
        idx_ref=[floor(n_echo_intl/2)+1+1:...
                 floor(n_echo_intl/2)+1+para.n_refs(1)];
        ref=dx(:,idx_ref(1)-1:idx_ref(end),:,:,:,:,:);
        ref_pol=ro_pol(idx_ref(1)-1:idx_ref(end));
        %
        disp('*** Striping reference echo information ***');
        dx(:,idx_ref,:,:,:,:,:)=[];
        ro_pol(idx_ref)=[];
        ro_pol_blpo(idx_ref)=[];
        te(idx_ref)=[];
        k(:,idx_ref,:)=[];
        disp('*** Take the average across partitions for ref data. This may not always work correctly! ***')
        ref=mean(ref,5);
    else
        ref=[];
        ref_pol=[];
    end
    
    
    if ~para.isgre && ...
            ~para.b_epi_positive && ~no_pc
        disp('*** EPI ghost correction ***');
        %% Need to process individual contrasts, mainly for multi-echo EPI
        if ~isempty(dblpo)
            for i=1:ncontr
                dx_contr=dx(:,(i-1)*n_echo_contr+1:i*n_echo_contr,:,:,:);
                if mid_blpo~=mid
                    dblpo_contr=dblpo;
                    ro_pol_blpo=ro_pol_blpo(1:n_echo_contr);
                else
                    dblpo_contr=dblpo(:,(i-1)*n_echo_contr+1:...
                                      i*n_echo_contr,:,:);
                    ro_pol_blpo=ro_pol((i-1)*n_echo_contr+1:i*n_echo_contr);
                end
                ro_pol_contr=ro_pol((i-1)*n_echo_contr+1:i*n_echo_contr);
                dx_contr=pha_crct_epi(dx_contr,dblpo_contr,...
                                      ro_pol_contr,ro_pol_blpo,...
                                      ord_pha_crct);
                dx(:,(i-1)*n_echo_contr+1:i*n_echo_contr,:,:,:)=dx_contr;
            end
        elseif ~isempty(ref)
            dx=pha_crct_epi(dx,ref,ro_pol,ref_pol,ord_pha_crct);
        end
    end
    % correct B0
    % get navigator data
    nav=[];
    if ~no_b0_crct
        disp('*** Global B0 correction ***');
        if para.b_nav_en
            [dnav,nav_te]=prep_nav(mid,1,para,no_nav_regress);
            nav.d=dnav;
            nav.te=nav_te;
        elseif para.b_ste_en
            if isempty(ste)
                if exist(fn_ste)
                    ste=read_data(fn_ste);
                else
                    ste=prep_ste(mid,'no_main',1);
                end
            end
            
        end
        dx=b0_crct(dx,nav,ste,te,para);
    end
    % correct fov
    if ~no_fov_crct
        disp('*** FOV correction ***');
        dx=fov_crct(dx,k,para);
    end
    % generate output
    % the new dimension is 
    % GRE: x,y,partition,slice,channel,necho,repetition
    % EPI: x,y,partition,slice,channel,necho or ncontr,repetition
    if para.isgre
        dx=reshape(dx,[nr,n_image_echo_shot,ns,nch,...
                       ninterl,npar/s2,n_delays*nreps]);
        dx=permute(dx,[1,5,6,3,4,2,7]);
        dx=reshape(dx,[nr,ninterl,...
                       npar/s2,ns,nch,...
                       n_image_echo_shot,n_delays*nreps]);
    else
        dx=reshape(dx,[nr,n_echo_contr,ncontr,ns,nch,...
                       ninterl,npar/s2,n_delays*nreps]);
        dx=permute(dx,[1,6,2,7,4,5,3,8]);
        dx=reshape(dx,[nr,ninterl*n_echo_contr,...
                       npar/s2,ns,nch,ncontr,n_delays*nreps]);
    end
    % apodization in phase-encoding directions
    dx=dx.*...
       reshape(win_tukey(np,apodiz,[0:s1:np-1]/np),...
               [1,np/s1]);
    if para.dimen==3
        dx=dx.*...
           reshape(win_tukey(npar,apodiz,...
                             [0:s2:npar-1]/npar),...
                   [1,1,npar/s2]);
    end
    % return k-space data
    if k_return
        % convert to k space
        y=fftmr(dx,1,1);
        return;
    end
    if no_comb
        % normalize channels based on noise measurement
        dx=covNorm(dx,cov_mat,5);
        y=fftmr(dx,-1,[2,3]);
    else
        if ~isempty(mid_b1)
            % Combine based on sensitivity profile
            % Work in progress:
            if exist(fn_b1)
                load(fn_b1);
            else
                b1info=gen_b1_info(mid_b1);
            end
            coord=get_coordinate(para,1,0);
            if ~isequal(size(b1info.coord),size(coord)) || ~isequalfp(coord,b1info.coord,1e-4)
                [sensit,smask]=interp_sense(b1info.coord,coord,b1info.sensit,b1info.mask);
            end
            % note that caipi is not supported for now.
            dk=fftmr(dx,1,1);
            y=zeros(nr,np,npar,ns);
            for is=1:ns
                y(:,:,:,is)=sense_recon(dk(:,:,:,is,:),sensit(:,:,is,:),smask(:,:,is),inv_cov,s1,s2,0,0,0);
            end
            disp('*** Performing SENSE reconstruction; CAIPI is not supported for now. ***');
            y=sense_recon(k,sensit,smask,inv_cov,s1,s2,0,0,0);
        else
            dx=covNorm(dx,cov_mat,5);
            y=fftmr(dx,-1,[2,3]);
% $$$             b1data=squeeze(y(:,:,:,1,:,1));
% $$$             m=auto_im_mask(sum(abs(b1data),4));
% $$$             [sensit,ref,~]=sense_m(b1data,m,eye(nch,nch));
% $$$             sensit=reshape(sensit,[nr,np,ns*npar,1,nch]);
% $$$             y=sum(y.*conj(sensit),5);
% $$$             y=abs(y).^0.5.*exp(1i*angle(y));
            y=mean(abs(y).^2,5).^0.5;            
    end

% $$$     end
end
