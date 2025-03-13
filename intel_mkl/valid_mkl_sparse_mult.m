function r=valid_mkl_sparse_mult(nch, n, m, verbose)
    if nargin<1
        nch=16;
    end
    if nargin<2
        n=1000;
    end
    if nargin<3
        m=200;
    end
    if nargin<4
        verbose=0;
    end
    if ischar(nch)
        nch=eval(nch);
    end
    if ischar(n)
        n=eval(n);
    end
    if ischar(m)
        m=eval(m);
    end
    if ischar(verbose)
        verbose=eval(verbose);
    end
    % x=randn(n,nch,'single')+1i*randn(n,nch,'single');
    x=squeeze(randn(n,nch,'double'));
    nnz_row=25;
    I=zeros(nnz_row*m,1);
    J=zeros(nnz_row*m,1);
    V=randn(nnz_row*m,1,'double');
    for i=1:m
        I((i-1)*nnz_row+1:i*nnz_row)=i;
        J((i-1)*nnz_row+1:i*nnz_row)=randi(n,nnz_row,1);
    end
    s=sparse(I,J,V);
    [val, row_ptr, col_ind] = sparse2csr(s, 0,'double');

    r=false;
    if verbose
        tic;
    end
    y=sparse_csr_mm(val,col_ind,row_ptr,m,n,x);
    if verbose
        toc;
    end
    yh=sparse_csr_mm_prit(val,col_ind,row_ptr,m,n,x.');


    %valid transpose
    n=3e2;
    m=n;
    nnz_row=1;
    I=[1:m];
    J=[1:m];
    V=ones(m,1)+0.5*1i;
    s=sparse(I,J,V,m,n);
    [v,r,c]=sparse2csr(s,0,'double');
    [vt,ct,rt]=mkl_sp_transpose_c16(double(v),c,r,int64(m),int64(n));
    r=true;

end