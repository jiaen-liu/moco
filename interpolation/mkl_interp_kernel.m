function ker=mkl_interp_kernel(coord_ori,coord_new,interp)
    if nargin<3
        interp='linear';
    end
    [irow,icol,val]=interp_kernel(coord_ori,coord_new,interp);
    n=numel(coord_new(:,:,:,1));
    base=0;
    col_ind=int64(icol-1+base);
    row_ptr=int64([0;cumsum(histc(irow,1:n))]+base);
    ker.col_ind=col_ind;
    ker.row_ptr=row_ptr;
    ker.val=val;
end