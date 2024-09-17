function [val, row_ptr, col_ind] = sparse2csr(A, base, data_type, varargin)
%SPARSE2CSR Convert MATLAB sparse matrix to CSR format.
%   To convert a CSC matrix instead of a CSR matrix to a sparse matrix
%   replace row_ptr with col_ptr, col_ind by row_ind, and transpose A. This
%   function supports both 0- an 1-based indexing
%
%   [val, row_ptr, col_ind] = SPARSE2CSR(A, base)
%
%   Inputs:
%   A        - MATLAB sparse matrix.
%  
%   base     - Either 0 or 1, for 0- or 1-based indexing of row_ptr and
%              col_ind.
%
%   Outputs:
%   val      - Vector of matrix entries.
%
%   row_ptr  - Pointer to index start of each row.
%
%   col_ind  - Index of elements in each colum.
%
%   References:
%   http://en.wikipedia.org/wiki/Sparse_matrix
%
%   Author:
%   Chris Fougner
    transpose=0;
    nvarargin=length(varargin);
    if nvarargin>0
        transpose=varargin{1};
    end
if nargin < 2 || isempty(base)
  base = 0;
end
if base ~= 0 && base ~= 1
  error('Base must be either 0 or 1.')
end


if transpose
    m = size(A, 2);
    [col, row, val] = find(A);
else
    m = size(A, 1);
    [col, row, val] = find(A.');
end

col_ind = int64(col - 1 + base);
row_ptr = int64([0; cumsum(histc(row, 1:m))] + base);
if strcmp(data_type, 'single')
    val=single(val);
end

end