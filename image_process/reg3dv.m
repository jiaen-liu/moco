function [im_reg,tra_mat,tra_par]=reg3dv(move,fix,varargin)
    p=inputParser;
    method='matlab';
    verbose=0;
    mcf_cost='normcorr';
    mask=[];
    fid=[];
    refine=[];
    nscale=0;
    inplane=0;
    res=[1,1,1];
    addParameter(p,'inplane',inplane,@isnumeric);
    addParameter(p,'method',method,@ischar);
    addParameter(p,'mcf_cost',mcf_cost,@ischar);
    addParameter(p,'res',res,@isnumeric);
    addParameter(p,'verbose',verbose,@isnumeric);
    addParameter(p,'mask',mask,@islogical);
    addParameter(p,'fid',fid,@isnumeric);
    addParameter(p,'refine',refine,@isnumeric);
    addParameter(p,'nscale',nscale,@isnumeric);
    p.parse(varargin{:});
    method=p.Results.method;
    mcf_cost=p.Results.mcf_cost;
    res=p.Results.res;
    verbose=p.Results.verbose;
    mask=p.Results.mask;
    fid=p.Results.fid;
    refine=p.Results.refine;
    nscale=p.Results.nscale;
    inplane=p.Results.inplane;
    si=size(move);
    ndim=numel(si);
    if ~isreal(move) || ~isreal(fix)
        absmove=abs(move);
        absfix=abs(fix);
    else
        absmove=move;
        absfix=fix;
    end
    ave_fix=mean(abs(fix(:)));
    absfix=absfix/ave_fix;
    absmove=absmove/ave_fix;
    nx=si(1);
    ny=si(2);
    if ndim>2
        nz=si(3);
    else
        nz=1;
    end
    
    if ndim == 3
        nv=1;
    elseif ndim==4
        nv=si(4);
    else
        nv=1;
        % error('The dimension should be 3 or 4');
    end

    
    im_reg=zeros(nx,ny,nz,nv);
    tra_mat=zeros(4,4,nv);
    tra_par=zeros(6,nv);
    [optimizer, metric] = imregconfig('monomodal');

    switch method
      case 'amri'
        fn_in=rp('im_coreg_amri_in.svd','shm');
        fn_out=rp('im_coreg_amri_out.svd','shm');
        fn_mot_par=rp('im_coreg_amri_mot_par.fit','shm');
        fn_mask=rp('im_coreg_amri_ref_mask.svd','shm');
        if ~isempty(fid)
            fid_char=int2str(fid);
            fn_in=file_addext(fn_in,fid_char);
            fn_out=file_addext(fn_out,fid_char);
            fn_mot_par=file_addext(fn_mot_par,fid_char);
            fn_mask=file_addext(fn_mask,fid_char);
        end;
        mag=single(cat(4,absmove,absfix));
        
        if exist(fn_in)
            delete(fn_in);
        end
        if exist(fn_out)
            delete(fn_out);
        end
        if exist(fn_mot_par)
            delete(fn_mot_par);
        end
        if exist(fn_mask)
            delete(fn_mask);
        end

        use_mask=~isempty(mask) && total(mask)>0;
        if use_mask
            mask=single(mask);
            mask(mask~=0)=1.0;
            save_data(fn_mask,mask);
        end

        save_data(fn_in,mag);
        cmd=['regist_svd ',...
             fn_out,' ',...
             fn_in,' ',...
             '-save ', fn_mot_par, ' '];
        if use_mask
            cmd=[cmd ' -mask ' fn_mask];
        end
        if inplane
            cmd=[cmd ' -inplane'];
        end
        system(cmd);
        tra_par=read_fitrec(fn_mot_par,[nx,ny,nz]);
        tra_par(:,end)=[];
        if nscale
            cmd=['regist_svd ',...
                 fn_out,' ',...
                 fn_in,' ',...
                 '-trans ', fn_mot_par, ' '];
        end
        im_reg=read_data(fn_out);
        im_reg(:,:,:,end)=[];
        im_reg=im_reg*ave_fix;
        for i=1:nv
            tra_mat(1:3,1:3,i)=rot3d(tra_par(1:3,i));
            tra_mat(1:3,4,i)=tra_par(4:6,i);
            tra_mat(4,4,i)=1;
        end
        if exist(fn_in)
            delete(fn_in);
        end
        if exist(fn_out)
            delete(fn_out);
        end
        if exist(fn_mot_par)
            delete(fn_mot_par);
        end
        if exist(fn_mask)
            delete(fn_mask);
        end
      case 'matlab'
        if verbose
            fprintf('The number of volumes to register: ');
        end
        for i=1:nv
            tmp=imregtform(absmove(:,:,:,i),absfix,'rigid',...
                           optimizer,metric);
            if isreal(move)
                im_reg(:,:,:,i)=imwarp(squeeze(move(:,:,:,i)),...
                                       imref3d([nx,ny,nz]),tmp,...
                                       'OutputView',imref3d([nx,ny,nz]),...
                                       'SmoothEdges',false);
            else
                im_reg(:,:,:,i)=imwarp(real(move(:,:,:,i)),imref3d([nx,ny,nz]),tmp,...
                                       'OutputView',imref3d([nx,ny,nz]),...
                                       'SmoothEdges',false)+...
                    1i*imwarp(imag(move(:,:,:,i)),imref3d([nx,ny,nz]),tmp,...
                              'OutputView',imref3d([nx,ny,nz]),...
                              'SmoothEdges',false);
            end
            tra_mat(:,:,i)=tmp.T;
            if verbose
                fprintf('%d ',nv-i);
            end
        end
        if verbose
            fprintf('\n');
        end
      case 'multires'
        if verbose
            fprintf('The number of volumes to register: ');
        end
        for i=1:nv
            tmp=estMotionMulti3(absfix,absmove(:,:,:,i));

            if isreal(move)
                im_reg(:,:,:,i)=warpAffine3(squeeze(move(:,:,:,i)),tmp,0);
            else
                im_reg(:,:,:,i)=warpAffine3(squeeze(real(move(:,:,:,i))),tmp,0)+...
                    1i*warpAffine3(squeeze(imag(move(:,:,:,i))),tmp,0);
            end
            tra_mat(:,:,i)=tmp;
            if verbose
                fprintf('%d ',nv-i);
            end
        end
        if verbose
            fprintf('\n');
        end
      case 'afni'
        fn_mov=rp('im_coreg_afni','shm');
        fn_out=rp('im_coreg_afni_out','shm');
        fn_fix=rp('im_coreg_afni_fix','shm');
        fn_mask=rp('im_coreg_afni_ref_mask.nii','shm');
        if ~isempty(fid)
            fid_char=int2str(fid);
            fn_mov=file_addext(fn_mov,fid_char);
            fn_out=file_addext(fn_out,fid_char);
            fn_fix=file_addext(fn_fix,fid_char);
            fn_mask=file_addext(fn_mask,fid_char);
        end;
        save_nii(make_nii(single(absmove)),[fn_mov '.nii']);
        save_nii(make_nii(single(absfix)),[fn_fix '.nii']);
        if exist([fn_out,'.nii'])
            delete([fn_out,'.nii']);
        end
        if exist(fn_mask)
            delete(fn_mask);
        end
        if exist([fn_mov,'.volreg_mats.aff12.1D'])
            delete([fn_mov,'.volreg_mats.aff12.1D']);
        end
        if exist([fn_mov,'.volreg_par'])
            delete([fn_mov,'.volreg_par']);
        end
        use_mask=~isempty(mask) && total(mask)>0;
        if use_mask
            mask=single(mask);
            mask(mask~=0)=1.0;
            save_nii(make_nii(mask),fn_mask);
        end
        if verbose
            str_verb=' -verbose ';
        else
            str_verb=' ';
        end
        if ~use_mask
        cmd=['3dvolreg', str_verb,...
             ' -1Dmatrix_save ',fn_mov,'.volreg_mats', ...
             ' -1Dfile ',fn_mov,'.volreg_par', ...
             ' -base ', fn_fix,'.nii''[0]''', ...
             ' -prefix ',fn_out,'.nii', ...
             ' ',fn_mov,'.nii'];
        else
            cmd=['3dvolreg', ' -weight ' fn_mask '''[0]''', str_verb,...
             ' -1Dmatrix_save ',fn_mov,'.volreg_mats', ...
             ' -1Dfile ',fn_mov,'.volreg_par', ...
             ' -base ', fn_fix,'.nii''[0]''', ...
             ' -prefix ',fn_out,'.nii', ...
             ' ',fn_mov,'.nii'];
        end
        cmd
        system(cmd);
        if exist([fn_out,'.nii'])
            im_tmp=load_nii([fn_out,'.nii']);
            im_reg=im_tmp.img;
        end
        if exist([fn_mov,'.volreg_mats.aff12.1D'])
            tra_mat=dlmread([fn_mov,'.volreg_mats.aff12.1D'],...
                            '',1,0);
        else
            tra_mat=repmat([1,0,0,0,...
                            0,1,0,0,...
                            0,0,1,0],...
                           [nv,1]);
        end

        tra_mat=tra_mat';
        tra_mat=reshape(tra_mat,[4,3,nv]);
        tra_mat=permute(tra_mat,[2,1,3]);
        tra_mat_tmp=zeros(4,4,nv);
        tra_mat_tmp(4,4,:)=1;
        tra_mat_tmp(1:3,:,:)=tra_mat;
        tra_mat=tra_mat_tmp;

        if nargout==3
            tra_par=dlmread([fn_mov,'.volreg_par'],...
                            '',0,0);
            tra_par=tra_par.';
            
        end
      case 'mcf'
        fn_mov=rp('im_coreg_mcf','shm');
        fn_out=rp('im_coreg_mcf_out','shm');
        fn_fix=rp('im_coreg_mcf_fix','shm');
        
        if ~isempty(fid)
            fid_char=int2str(fid);
            fn_mov=file_addext(fn_mov,fid_char);
            fn_out=file_addext(fn_out,fid_char);
            fn_fix=file_addext(fn_fix,fid_char);
        end;
        
        if exist([fn_out,'.nii.gz'])
            delete([fn_out,'.nii.gz']);
        end
        if exist([fn_mov,'.nii'])
            delete([fn_mov,'.nii']);
        end
        if exist([fn_fix,'.nii'])
            delete([fn_fix,'.nii']);
        end
        for i=1:nv
            tmp=[fn_out '.nii.gz.mat' filesep 'MAT_' sprintf('%04d',i-1)];
            if exist(tmp)
                delete(tmp);
            end
        end
        if exist([fn_out,'.nii.gz.mat'])
            rmdir([fn_out,'.nii.gz.mat'],'s');
        end        
        niitmp=make_nii(absmove);
        niitmp.hdr.dime.pixdim(1)=length(res);
        niitmp.hdr.dime.pixdim(2:1+length(res))=res;
        save_nii(niitmp,[fn_mov '.nii']);
        niitmp=make_nii(absfix);
        niitmp.hdr.dime.pixdim(1)=length(res);
        niitmp.hdr.dime.pixdim(2:1+length(res))=res;
        save_nii(niitmp,[fn_fix '.nii']);
        if verbose
            str_verb=' -report';
        else
            str_verb=' ';
        end
        
        cmd_mcf=['mcflirt -in ' fn_mov '.nii' ...
                 ' -o ' fn_out, '.nii.gz' ...
                 ' -r ',fn_fix '.nii' ...
                 ' -cost ' mcf_cost ' -sinc_final'...
                 ' -mats' str_verb];
        disp(cmd_mcf);
        system(cmd_mcf);
        im_tmp=load_nii([fn_out '.nii.gz']);
        im_reg=im_tmp.img;
        im_reg=im_reg*ave_fix;
        for i=1:nv
            tra_mat(:,:,i)=dlmread([fn_out '.nii.gz.mat' filesep 'MAT_' sprintf('%04d',i-1)]);
        end
    end
    if ~isempty(refine) && refine==1
        std_coreg=std(im_reg,0,4);
        if isempty(mask)
            mask_tmp=true(nx,ny,nz);
        else
            mask_tmp=mask;
        end
        std_thrd=prctile(std_coreg(logical(mask_tmp)),75);
        mask_std=std_coreg<std_thrd;
        if total(mask_std)==0
            mask_std=mask_tmp;
        end
        mask_tmp=mask_tmp&mask_std;

        [im_reg,tra_mat,tra_par]=reg3dv(move,fix,'method',method,...
                                        'mask',mask_tmp,...
                                        'fid',fid);
    end
        
end