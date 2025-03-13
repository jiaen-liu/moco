function y=sparse_csr_mm(val,col_ind,row_ptr,m,n,x)
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
            y=sparseMultiMatC8(val,col_ind,row_ptr,m,n,x);
        elseif ~isreal(val) && isreal(x)
            y=sparseMultiMatSingle(real(val),col_ind,row_ptr,m,n,x)+...
              1i*sparseMultiMatSingle(imag(val),col_ind,row_ptr,m,n,x);
        elseif isreal(val) && ~isreal(x)
            y=sparseMultiMatSingle(val,col_ind,row_ptr,m,n,real(x))+...
              1i*sparseMultiMatSingle(val,col_ind,row_ptr,m,n,imag(x));
        else
            % both real
            y=sparseMultiMatSingle(val,col_ind,row_ptr,m,n,x);
        end
    elseif isa(val,'double') && isa(x,'double')
        % use double functions
        if ~isreal(val) && ~isreal(x)
            y=sparseMultiMatC16(val,col_ind,row_ptr,m,n,x);
        elseif ~isreal(val) && isreal(x)
            y=sparseMultiMatDouble(real(val),col_ind,row_ptr,m,n,x)+...
              1i*sparseMultiMatDouble(imag(val),col_ind,row_ptr,m,n,x);
        elseif isreal(val) && ~isreal(x)
            y=sparseMultiMatDouble(val,col_ind,row_ptr,m,n,real(x))+...
              1i*sparseMultiMatDouble(val,col_ind,row_ptr,m,n,imag(x));
        else
            % both real
            y=sparseMultiMatDouble(val,col_ind,row_ptr,m,n,x);
        end
    else
        val=single(val);
        x=single(x);
        % use single functions
        if ~isreal(val) && ~isreal(x)
            y=sparseMultiMatC8(val,col_ind,row_ptr,m,n,x);
        elseif ~isreal(val) && isreal(x)
            y=sparseMultiMatSingle(real(val),col_ind,row_ptr,m,n,x)+...
              1i*sparseMultiMatSingle(imag(val),col_ind,row_ptr,m,n,x);
        elseif isreal(val) && ~isreal(x)
            y=sparseMultiMatSingle(val,col_ind,row_ptr,m,n,real(x))+...
              1i*sparseMultiMatSingle(val,col_ind,row_ptr,m,n,imag(x));
        else
            % both real
            y=sparseMultiMatSingle(val,col_ind,row_ptr,m,n,x);
        end
    end
end