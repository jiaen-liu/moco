function [c,rp]=coo2csr(ic,ir,m,base,nnz_row)
% the input should be ordered in row major mode
c = int64(ic - 1 + base);
if nargin>4
    % a fast way if know the nnz per row
    rp=int64([0:m].'*nnz_row+base);
else
    rp = int64([0; cumsum(histc(ir, 1:m))] + base);
end
end