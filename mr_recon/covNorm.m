function d=covNorm(d,cov_mat,dim)
% last dimention of d is channel
    si=size(d);
    ndim=length(si);
    if nargin<3
        dim=ndim;
    end
    idx_trps=[1:ndim];
    idx_trps(end)=dim;
    idx_trps(dim)=ndim;
    if dim<ndim
        d=permute(d,idx_trps);
    end
    si_trps=si;
    si_trps(end)=si(dim);
    si_trps(dim)=si(end);
    nch=si(dim);
    n=numel(d)/nch;
    d=reshape(d,[n,nch]);
    inv_cov=inv(cov_mat);
    chol_mat=conj(chol(inv_cov,'lower'));
    d=d*chol_mat;
    d=reshape(d,si_trps);
    if dim<ndim
        d=permute(d,idx_trps);
        d=reshape(d,si);
    end
end