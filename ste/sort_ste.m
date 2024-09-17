% 2022-12-23: Jiaen Liu, added support for Philips
function [y]=sort_ste(data)
    if ~isfield(data,'vendor')
        vendor='siemens';
    else
        vendor=data.vendor;
        if iscell(vendor)
            vendor=vendor{1};
        end
    end
  data=struct2double(data);
  % icontr=para.icontr;
  ncontr=data.pe.ncontr;
  % parellel imaging parameters
  s1=data.pe.sp;
  s2=data.pe.r/s1;
  [nr,necho,nch,nshot]=size(data.d);
  np=data.para.steref_dim_p;
  ns=data.para.steref_dim_s;
  n=nr*np*ns;
  n_nch=n*nch;
  % inverse covariance matrix
  inv_cov=inv(data.cov_mat);
  cov_mat=(data.cov_mat);
  % number of echoes per contrast
  necho_contr=necho/ncontr;
  % r is the number of partial volumes to cycle through a full
  % fov; r may depend on the ste parameters.
  % r=s1^double(data.pe.y_cyc>0)*s2^data.pe.z_cyc;
  r=s1*s2;
  dkz=data.pe.dkz;
  pe_dir=ste_pe_dir(data);
  % 
  nshot_v=ns/s2;
  nshot_fv=nshot_v*r;
  % find the blipless shots
  k=reshape(data.pe.k,[2*necho,data.pe.hl(1)/necho]);
  mblpl=(sum(abs(k),1) == 0);
  k=reshape(data.pe.k,[2,necho,data.pe.hl(1)/necho]);
  % non-blipless k coordinates
  k_nblpl=k(:,:,~mblpl);
  te=double(data.te);
  % blipless interval
  dblpl=data.pe.hl(1)/necho/sum(mblpl);
  % non-blipless data
  d_nblpl=data.d;
  mblpl=mblpl(:);
  nblpl_fv=total(mblpl);
  mblpl=repmat(mblpl,[ceil(nshot/(data.pe.hl(1)/necho)),1]);
  mblpl=mblpl(1:nshot);
  
  d_nblpl(:,:,:,mblpl)=[];
  % non-blipless shots
  nshot_nblpl=size(d_nblpl,4);
  % acquired partial volumnes
  nv=floor(nshot_nblpl/nshot_v);
  nvf=floor(nshot_nblpl/nshot_v/r);
  % shot index
  % Each colume contains:
  % index of full fov, partial fov and blipless ste
  idxshot=zeros(3,nshot);
  iblpl=0;
  iparfov=0;
  ifulfov=0;
  % mark the index of non-blipless shots
  nshot_nblpl=sum(~mblpl);
  idx_nblpl=zeros(nshot,1);
  idx_nblpl(~mblpl)=[1:nshot_nblpl];
  % mark the index of shots where the center of k-space is sampled
  idx_shot_cenff=zeros(nvf,1);
  idx_shot_cenpf=zeros(nv,1);
  for i=1:nshot
      itmp=floor((idx_nblpl(i)-1)/nshot_fv)+1;
      if mod(idx_nblpl(i)-1, nshot_fv)== floor(nshot_fv/2)...
              && itmp<=nvf
          idx_shot_cenff(itmp)=i;
      end
      itmp=floor((idx_nblpl(i)-1)/nshot_v)+1;
      if mod(idx_nblpl(i)-1, nshot_v)== floor(nshot_v/2)...
              && itmp<=nv
          idx_shot_cenpf(itmp)=i;
      end
      if ~mblpl(i)
          iparfov=iparfov+1;
          ifulfov=ifulfov+1;
      else
          iblpl=iblpl+1;
      end
      % mark the index of navigator when each shot is taken
      idx_shot(1,i)=min(max(0,floor((ifulfov-1)/nshot_fv)+1),nvf);
      idx_shot(2,i)=min(max(0,floor((iparfov-1)/nshot_v)+1),nv);
      idx_shot(3,i)=iblpl;
      % blipless scan belongs to the comming volume
      if mblpl(i) && mod(ifulfov-1,nshot_fv) == (nshot_fv-1)
          idx_shot(1,i)=min(idx_shot(1,i)+1,nvf);
      end
      if mblpl(i) && mod(iparfov-1,nshot_v) == (nshot_v-1)
          idx_shot(2,i)=min(idx_shot(2,i)+1,nv);
      end
  end
  y=struct('d_nblpl',d_nblpl,...
           'idx_nblpl',idx_nblpl,...
           'idx_shot_cenff',idx_shot_cenff,...
           'idx_shot_cenpf',idx_shot_cenpf,...
           'idx_shot',idx_shot,...
           'nshot_fv',nshot_fv,...
           'nshot_v',nshot_v,...
           'ncontr',ncontr,...
           'k',k,...
           'k_nblpl',k_nblpl,...
           'dims_k',[nr,np,ns,nch],...
           'y_cyc',data.pe.y_cyc,...
           'z_cyc',data.pe.z_cyc,...
           'nparv_cyc',r,...
           'inv_cov',inv_cov,...
           'para',data.para,...
           'mblpl',mblpl,...
           'df_blpl',data.df_blpl(:,2),...
           'vendor',vendor);
  if ~isfield(y.para,'nk_shot')
      y.para.nk_shot=1;
  end
  if ~isfield(y.para,'dkz_caipi')
      y.para.dkz_caipi=0;
  end
  if isfield(data,'k_pimg') && ~isempty(data.k_pimg)
      y=setfield(y,'k_pimg',data.k_pimg);
      y=setfield(y,'para_pimg',data.para_pimg);
      y=setfield(y,'mid_pimg',data.mid_pimg);
  end
  if isfield(data,'sense_philips')
      y=setfield(y,'sense_philips',data.sense_philips);
  end
  if field_true(data,'combined')
      y=setfield(y,'combined',1);
      y=setfield(y,'ishot_break',data.ishot_break);
  else
      y=setfield(y,'combined',0);
  end
  if isfield(data,'kyz')
      y=setfield(y,'kyz',data.kyz);
  end
  ytmp=var2struct('s1','s2','necho','nshot','dkz','nshot_nblpl',...
                  'nv','nvf','te','nblpl_fv','pe_dir','cov_mat');
  fldn_tmp=fieldnames(ytmp);
  for i=1:numel(fldn_tmp)
      y=setfield(y,fldn_tmp{i},getfield(ytmp,fldn_tmp{i}));
  end
end
