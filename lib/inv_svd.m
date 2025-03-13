function inv=inv_svd(A)
    [u,s,v]=svd(A);
    s_diag=s(idx_diag(s));
    if s_diag(1)~=0
        idx=find(s_diag>s_diag(1)*1e-3);
        idxz=find(s_diag<=s_diag(1)*1e-3);
        s_diag(idx)=1./s_diag(idx);
        s_diag(idxz)=0;
    end
    s2=s;
    s2(idx_diag(s))=s_diag;
    inv=(u*s2*v')';
end