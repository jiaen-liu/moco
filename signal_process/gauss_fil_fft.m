function y=gauss_fil_fft(data, tao, dt, dim)
    data_cp = data;
    si = size(data_cp);
    if nargin<4
        dim=find(si~=1,1,'first');
    end


    trp_idx = [1:numel(si)];
    trp_idx(1) = dim;
    trp_idx(dim) = 1;
    
    
    data_cp = permute(data_cp, trp_idx);
    si_trp = size(data_cp);
    
    sig = log(2)^0.5/1.414/pi*tao/dt;
    width = round(6*sig);
    if mod(width,2) == 0
        width = width+1;
    end
    hwidth = floor(width/2);
    gauss = normpdf([-hwidth:hwidth], 0, sig);
    gauss=gauss/total(gauss);

    l = si_trp(1);
    ll = 2^nextpow2(l+width-1);
    if numel(si_trp) > 1
        tmp = zeros([ll, si_trp(2:end)]);
    else
        tmp = zeros(ll, 1);
    end
    
    tmp(ll/2-floor(l/2)+1:ll/2+floor((l-1)/2)+1, :,:,:,:,:,:,:) = data_cp;
    kernel = zeros(ll,1);
    kernel(ll/2-hwidth:ll/2+hwidth) = gauss;

    if numel(si_trp) > 1
        kernel = repmat(kernel, [1, si_trp(2:end)]);
    end
    
    tmp = (fftmr(fftmr(kernel, -1, 1).*fftmr(tmp, -1, 1), 1, 1))*ll;
    y=permute(tmp(ll/2-floor(l/2):ll/2+floor((l-1)/2), :,:,:,:,:,:,:), trp_idx);
end
