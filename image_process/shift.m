function y=shift(a,k)
    si=size(a);
    ndim=numel(si);
    ndim_k=numel(size(k));
    if ndim_k>ndim
        error('The dimension of the array should be no more than that of the shift parameter!');
    end
    idx_nnz=find(k);

    y=a;
    for i=1:numel(idx_nnz)
        ktmp=k(idx_nnz(i));
        if rem(ktmp,1)==0
            y=circshift(y,ktmp,idx_nnz(i));
        else
            n=si(idx_nnz(i));
            y=fftmr(y,-1,idx_nnz(i));
            phase_ramp=exp(-1i*2*pi*ktmp*((0:n-1)/n-0.5));
            if ~mod(n,2)
                phase_ramp(1)=real(phase_ramp(1));
            end
            si_reshape=ones(1,ndim);
            si_reshape(idx_nnz(i))=n;
            y=y.*reshape(phase_ramp,si_reshape);
            y=fftmr(y,1,idx_nnz(i));
        end
    end
    if isreal(a)
        y=real(y);
    end
end