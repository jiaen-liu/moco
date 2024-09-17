function y=idx_diag(A)
    [m,n]=size(A);
    y=[1:min(m,n)]*(m+1)-m;
    y=y(:);
end
