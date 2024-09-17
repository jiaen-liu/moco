function y=regridding_arr(d,para,apodization)
    si=size(d);
    d=reshape(d,[si(1),prod(si(2:end))]);
    ftmat=calc_ft_mat(para,apodization);
    y=ftmat*d;
    si(1)=size(y,1);
    y=reshape(y,si);
end