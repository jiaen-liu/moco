function y=sos(x)
    dims=size(x);
    y=sum(abs(x).^2,numel(dims)).^0.5;
end