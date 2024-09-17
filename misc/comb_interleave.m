function c=comb_interleave(a,b)
    a=reshape(a,[1,numel(a)]);
    b=reshape(b,[1,numel(b)]);
    if numel(a)>numel(b)
        c=[a(1:end-1);b];
        c=[c(:);a(end)];
    else
        c=[a;b];
        c=c(:);
    end
end