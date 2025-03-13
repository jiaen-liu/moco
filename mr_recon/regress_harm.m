function do=regress_harm(d,n,period,order,period_sample)
% period: period of the source
% order: number of harmonics to remove
% period_sample: period of the sample (d)
% n: time points of the source
    mask=false(n,1);
    mask(1:period_sample:end)=true;
    % determine order if the matrix is low rank
    while true
        disp(num2str(order,'*** Try navigator regression with %dth order ***'));
        x=mtxHarmRegr(n,1,period,order,mask);
        % if rank(x)<min(size(x))
        if cond(x)>10
            order=order-1;
        else
            break;
        end
        if order==1
            warning('*** Only 1st order used in navigator regression! ***');
            break;
        end
    end
    
    c=multiRegress(d,x);
    do=d-x(:,2:end)*c(2:end,:);
end
