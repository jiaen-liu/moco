function do=move_matrix(d,mot,direction,varargin)
    si=size(d);
    if numel(si)==4
        nv=si(4);
        nz=si(3);
    elseif numel(si)==2
        nv=1;
        nz=1;
    elseif numel(si)==3
        nv=1;
        nz=si(3);
    end
    nx=si(1);
    ny=si(2);
    fid=[];
    n=nx*ny*nz;
    method='matlab';
    interp='cubic';
    p=inputParser;
    res=[1,1,1];
    addParameter(p,'interp',interp,@ischar);
    addParameter(p,'res',res,@isnumeric);
    addParameter(p,'method',method,@ischar);
    addParameter(p,'fid',fid,@isnumeric);
    p.parse(varargin{:});
    interp=p.Results.interp;
    res=p.Results.res;
    res=res(:).';
    method=p.Results.method;
    fid=p.Results.fid;
    mot=mot(:);
    breal=isreal(d);
    if breal
        do=zeros(nx,ny,nz,nv);
    else
        do=complex(zeros(nx,ny,nz,nv));
    end
    switch method
      case 'matlab'
        R=rot3d(mot(1:3));
        % the matlab x and y axises are switched.
        Rtmp=R;
        R(:,1)=Rtmp(:,2);
        R(:,2)=Rtmp(:,1);
        Rtmp=R;
        R(1,:)=Rtmp(2,:);
        R(2,:)=Rtmp(1,:);
        T=mot([5,4,6]);
        tra_matrix=[[R,T]',[0;0;0;1]];

        if strcmp(direction,'inverse')
            tra_matrix=inv(tra_matrix);
        end
        tra_matrix(:,4)=[0;0;0;1];
        
        if nz==1
            imref=@imref2d;
            tra_matrix2d=zeros(3,3);
            tra_matrix2d(1:2,1:2)=tra_matrix(1:2,1:2);
            tra_matrix2d(3,1:2)=tra_matrix(4,1:2);
            tra_matrix2d(3,3)=1;
            tra_matrix=tra_matrix2d;
            tform=affine2d(eye(3,3));
            Rin=imref([nx,ny,nz],res(2),res(1));
            Rin.XWorldLimits=Rin.XWorldLimits-mean(Rin.XWorldLimits);
            Rin.YWorldLimits=Rin.YWorldLimits-mean(Rin.YWorldLimits);
        else
            imref=@imref3d;
            tform=affine3d(eye(4,4));
            Rin=imref([nx,ny,nz],res(2),res(1),res(3));
            Rin.XWorldLimits=Rin.XWorldLimits-mean(Rin.XWorldLimits);
            Rin.YWorldLimits=Rin.YWorldLimits-mean(Rin.YWorldLimits);
            Rin.ZWorldLimits=Rin.ZWorldLimits-mean(Rin.ZWorldLimits);
        end
        Rout=Rin;
        tform.T=tra_matrix;
        for i=1:nv
            if breal
                do(:,:,:,i)=imwarp(d(:,:,:,i),...
                                   Rin,...
                                   tform,...
                                   interp,...
                                   'OutputView',Rout);
            else
                do(:,:,:,i)=imwarp(real(d(:,:,:,i)),...
                                   Rin,...
                                   tform,...
                                   interp,...
                                   'OutputView',Rout)+...
                    1i*imwarp(imag(d(:,:,:,i)),...
                              Rin,...
                              tform,...
                              interp,...
                              'OutputView',Rout);
            end
        end
      case 'flirt'
        % convert motion par to affine matrix in fsl format
        if nv>1
            error('*** This program based on flirt only transforms single volume at a time! ***');
        end
        R=rot3d(mot(1:3));
        T=a2v(mot([4,5,6]));
        tra_matrix=[[R,T];[0,0,0,1]];
        tra_matrix=convert_afn(tra_matrix,...
                               eye(3,3),...
                               [-1,1,1],...
                               -[-nx+1,ny-1,nz-1].*res/2);
        if strcmp(direction,'inverse')
            tra_matrix=inv(tra_matrix);
        end
        tra_matrix(4,:)=[0,0,0,1];
        fn_pre=fullfile('~','result','shm','move_matrix_flirt');
        if ~isempty(fid)
            fn_pre=file_addext(fn_pre,['_fid', num2str(fid)]);
        end
        fn_in=[fn_pre,'_in.nii.gz'];
        fn_ref=[fn_pre,'_ref.nii.gz'];
        fn_out=[fn_pre,'_out.nii.gz'];
        fn_mat_folder=[fn_pre,'.mat'];
        fn_mat=fullfile(fn_mat_folder,'MAT_00000');
        
        cmd=['flirt -in ' fn_in ' -ref ' fn_ref ...
             ' -out ' fn_out ' -init ' fn_mat ...
             ' -applyxfm -interp spline'];
        % disp(cmd);
        if ~breal
            dr=real(d);
        else
            dr=d;
        end
        % transform the real part
        niitmp=make_nii(dr);
        niitmp.hdr.dime.pixdim(1)=length(res);
        niitmp.hdr.dime.pixdim(2:1+length(res))=res;
        save_nii(niitmp,fn_in);
        
        niitmp=make_nii(dr(:,:,:,1));
        niitmp.hdr.dime.pixdim(1)=length(res);
        niitmp.hdr.dime.pixdim(2:1+length(res))=res;
        save_nii(niitmp,fn_ref);
        save_afnmat_fsl(fn_mat_folder,tra_matrix);
        status=system(cmd);
        if status~=0
            error('*** Error transforming matrix using flirt! ***');
        end
        niitmp=load_nii(fn_out);
        do=niitmp.img;
        if ~breal
            % transform the imagery part
            di=imag(d);
            niitmp=make_nii(di);
            niitmp.hdr.dime.pixdim(1)=length(res);
            niitmp.hdr.dime.pixdim(2:1+length(res))=res;
            save_nii(niitmp,fn_in);
            
            niitmp=make_nii(di(:,:,:,1));
            niitmp.hdr.dime.pixdim(1)=length(res);
            niitmp.hdr.dime.pixdim(2:1+length(res))=res;
            save_nii(niitmp,fn_ref);
            status=system(cmd);
            if status~=0
                error('*** Error transforming matrix using flirt! ***');
            end
            niitmp=load_nii(fn_out);
            do=do+1i*niitmp.img;
        end
        % clear the temporal files
        delete(fn_in);
        delete(fn_ref);
        delete(fn_out);
        rmdir(fn_mat_folder,'s');
    end
        

end