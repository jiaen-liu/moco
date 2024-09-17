function y=sparse_csr_mm_prit(val,col_ind,row_ptr,m,n,x)
    m=int64(m);
    n=int64(n);
    if ~isa(col_ind,'int64')
        col_ind=int64(col_ind);
    end
    if ~isa(row_ptr,'int64')
        row_ptr=int64(row_ptr);
    end
    
    if isa(val,'single') && isa(x,'single')
        % use single functions
        if ~isreal(val) && ~isreal(x)
            y=sparseMultiMatC8PriT(val,col_ind,row_ptr,m,n,x);
        elseif ~isreal(val) && isreal(x)
            y=sparseMultiMatSinglePriT(real(val),col_ind,row_ptr,m,n,x)+...
              1i*sparseMultiMatSinglePriT(imag(val),col_ind,row_ptr,m,n,x);
        elseif isreal(val) && ~isreal(x)
            y=sparseMultiMatSinglePriT(val,col_ind,row_ptr,m,n,real(x))+...
              1i*sparseMultiMatSinglePriT(val,col_ind,row_ptr,m,n,imag(x));
        else
            % both real
            y=sparseMultiMatSinglePriT(val,col_ind,row_ptr,m,n,x);
        end
    elseif isa(val,'double') && isa(x,'double')
        % use double functions
        if ~isreal(val) && ~isreal(x)
            y=sparseMultiMatC16PriT(val,col_ind,row_ptr,m,n,x);
        elseif ~isreal(val) && isreal(x)
            y=sparseMultiMatDoublePriT(real(val),col_ind,row_ptr,m,n,x)+...
              1i*sparseMultiMatDoublePriT(imag(val),col_ind,row_ptr,m,n,x);
        elseif isreal(val) && ~isreal(x)
            y=sparseMultiMatDoublePriT(val,col_ind,row_ptr,m,n,real(x))+...
              1i*sparseMultiMatDoublePriT(val,col_ind,row_ptr,m,n,imag(x));
        else
            % both real
            y=sparseMultiMatDoublePriT(val,col_ind,row_ptr,m,n,x);
        end
    else
        val=single(val);
        x=single(x);
        % use single functions
        if ~isreal(val) && ~isreal(x)
            y=sparseMultiMatC8PriT(val,col_ind,row_ptr,m,n,x);
        elseif ~isreal(val) && isreal(x)
            y=sparseMultiMatSinglePriT(real(val),col_ind,row_ptr,m,n,x)+...
              1i*sparseMultiMatSinglePriT(imag(val),col_ind,row_ptr,m,n,x);
        elseif isreal(val) && ~isreal(x)
            y=sparseMultiMatSinglePriT(val,col_ind,row_ptr,m,n,real(x))+...
              1i*sparseMultiMatSinglePriT(val,col_ind,row_ptr,m,n,imag(x));
        else
            % both real
            y=sparseMultiMatSinglePriT(val,col_ind,row_ptr,m,n,x);
        end
    end
end