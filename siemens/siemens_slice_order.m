function idx=siemens_slice_order(n)
idx1=[1:floor(n/2)];
idx2=[floor(n/2)+1:n];
if mod(n,2)==1
    % odd
    idx=[col([idx2(1:end-1);idx1]);idx2(end)];
else
    % even
    idx=col([idx2;idx1]);
end
end