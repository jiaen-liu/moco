function [vt,ct,rt]=mex_sp_transpose(v,c,r,m,n);
    if m~=length(r)-1
        error('*** The number of the row pointer elements should be equal to the number of rows plus one! ***');
    end
    if ~isa(c,'int64')
        c=int64(c);
    end
    if ~isa(r,'int64')
        r=int64(r);
    end
    
    if isa(v,'double')
        if isreal(v)
            % double
            [vt,ct,rt]=sp_transpose_double(v,c,r,int64(m),int64(n));
        else
            % double complex
            [vt,ct,rt]=sp_transpose_c16(v,c,r,int64(m),int64(n));
        end
    elseif isa(v,'single')
        if isreal(v)
            % single
            [vt,ct,rt]=sp_transpose_single(v,c,r,int64(m),int64(n));
        else
            % single complex
            [vt,ct,rt]=sp_transpose_c8(v,c,r,int64(m),int64(n));
        end
    end
end