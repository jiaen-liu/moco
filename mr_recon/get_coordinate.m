function p=get_coordinate(para,oversamp,steref,nonsiemens)
    if nargin<2
        oversamp=1;
    end
    if nargin<3
        steref=0;
    end
    if nargin==4 && nonsiemens
        nonsiemens=0;
        error('*** Non-Siemens coordinates are not supported ***');
    end
    if steref==1
        oversamp=1;
    end
    
    if ~field_true(para,'b_ste_en')
        steref=0;
    end
    
    if field_true(para,'b_ste_en')
        nr=conditional(steref,para.steref_dim_r,para.nr);
        np=conditional(steref,para.steref_dim_p,para.np);
        npar=conditional(steref,para.steref_dim_s,...
                         conditional(oversamp,para.n_partitions,...
                                     para.n_partitions_nos));
        resr=conditional(steref,para.steref_res_r,para.resr);
        resp=conditional(steref,para.steref_res_p,para.resp);
        ress=conditional(steref,para.steref_res_s,para.ress);
    else
        nr=para.nr;
        np=para.np;
        npar=conditional(oversamp,para.n_partitions,...
                         para.n_partitions_nos);
        resr=para.resr;
        resp=para.resp;
        ress=para.ress;
    end
    
    ns=para.n_slices*npar;
    p=zeros(nr,np,ns,3);

    pr=([0:nr-1]-(nr-1)/2)*resr;
    pp=([0:np-1]-(np-1)/2)*resp;
    part_shift=([0:npar-1]-(npar-1)/2)*ress;
    for is=1:para.n_slices
        new_coor=gen_coordinate(para.snormal(:, is), para.prot);
        for ipart=1:npar
            for ip=1:np
                for ir=1:nr
                    p(ir,ip,(is-1)*npar+ipart,:)=...
                        pr(ir)*new_coor(:,1)+...
                        pp(ip)*new_coor(:,2)+...
                        part_shift(ipart)*new_coor(:,3);
                end
            end
        end
        p(:,:,(is-1)*npar+1:is*npar,1)=...
            p(:,:,(is-1)*npar+1:is*npar,1)+para.x_shift(is);
        p(:,:,(is-1)*npar+1:is*npar,2)=...
            p(:,:,(is-1)*npar+1:is*npar,2)+para.y_shift(is);
        p(:,:,(is-1)*npar+1:is*npar,3)=...
            p(:,:,(is-1)*npar+1:is*npar,3)+para.z_shift(is);
    end
end