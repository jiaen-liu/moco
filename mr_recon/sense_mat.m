function [m,varn]=sense_mat(s, inv_cov, sp, ss, ds, dy, dz, slice_noise)
    s=double(s);
    inv_cov=double(inv_cov);
    [nr, np, ns, nch] = size(s);
    if nargin<8
        slice_noise=[1:ns];
    end
    n=nr*np*ns;
    [mcai,I,J,V]=caipi_mat(nr,np,ns,sp,ss,ds,dy,dz);
    % inverse of the noise covariance matrix
    sn=reshape(s,[n,nch])*conj(inv_cov);
    sh=reshape(conj(s),[n,nch]);
    m1=sparse([1:n],[1:n],sn(:,1),n,n);
    m2=sparse([1:n],[1:n],sh(:,1),n,n);
    m=m2*mcai*m1;
    if nch > 1 
        for ich=2:nch
            m1=sparse([1:n],[1:n],sn(:,ich),n,n);
            m2=sparse([1:n],[1:n],sh(:,ich),n,n);
            m=m+m2*mcai*m1; 
        end
    end
    
    if nargout>1
        % calculate noise
        varn=zeros(n,1);
        r=sp*ss;
        nd=n/r;
        mask=false(n,1);
        ivox=[1:nr*np].'+(nr*np)*(slice_noise(:).'-1);
        ivox=ivox(:);
        for i=1:length(ivox)
            if ~mask(ivox(i))
                ml=full(m(J(:,ivox(i)),J(:,ivox(i))));
                ml_inv=inv_svd(ml);
                test=ml_inv(idx_diag(ml_inv));
                if any(abs(test(:))>1e25)
                    a=1;
                end
                varn(J(:,ivox(i)),1)=ml_inv(idx_diag(ml_inv));
                mask(J(:,ivox(i)))=true;
            end
        end
        varn=reshape(varn,[nr,np,ns]);
    end
    
end