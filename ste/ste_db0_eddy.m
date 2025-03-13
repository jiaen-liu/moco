function db0=ste_db0_eddy(para,m1,fn_calib)
    load(fn_calib);
    % gradient of calibrated field difference
    % in the coordinate of the magnetic
    % for siemenst it is
    % ->x
    
    %  | y
    %  \/
    
    %  ^ z
    % /

    if numel(c)==3
        g=c;
    elseif numel(c)==4;
        g=c(2:4);
    end
    % transformation due to motion: m1
    % transformation from magnet frame to head frame
    [m0,~]=gen_coordinate(para.snormal,para.prot);
    t0=[para.x_shift;para.y_shift;para.z_shift];

    res=diag([para.steref_res_r,para.steref_res_p,para.steref_res_s]);



    nv=size(m1,3);
    db0=zeros(para.steref_dim_r,para.steref_dim_p,para.steref_dim_s,nv);
    % db02=zeros(para.steref_dim_r,para.steref_dim_p,para.steref_dim_s,nv);
    nx=para.steref_dim_r;
    ny=para.steref_dim_p;
    nz=para.steref_dim_s;

    x=[0:nx-1]-(nx)/2;
    y=[0:ny-1]-(ny)/2;
    z=[0:nz-1]-(nz)/2;

    [X,Y,Z]=ndgrid(x,y,z);

    r2=[X(:).';Y(:).';Z(:).'];
% $$$     tformtmp=affine3d(eye(4,4));
% $$$     coord=get_coordinate(mid,1,1);
% $$$     db02_tmp=g(1)*coord(:,:,:,1)+...
% $$$              g(2)*coord(:,:,:,2)+...
% $$$              g(3)*coord(:,:,:,3);
    for i=1:nv
        db0(:,:,:,i)=reshape(g.'*(m0*res*(m1(1:3,1:3,i)*r2+m1(1:3,4,i))+t0),...
                              [nx,ny,nz]);


        
% $$$         tformtmp.T=m1(:,:,i);
% $$$         
% $$$         db02(:,:,:,i)=imwarp(db02_tmp,...
% $$$                              imref3d([nx,ny,nz]),...
% $$$                              tformtmp,...
% $$$                              'OutputView',imref3d([nx,ny,nz]));
    end
end
