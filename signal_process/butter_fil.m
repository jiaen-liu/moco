function y=butter_fil(data, fc, fs, dim, order, ftype)
    si = size(data);
    if nargin<4
        dim=find(si~=1,1,'first');
    end
    if nargin<5
        order=6;
    end
    if nargin<6
        ftype=conditional(numel(fc)==2,'bandpass','low');
    end
    
    trp_idx = [1:numel(si)];
    trp_idx(1) = dim;
    trp_idx(dim) = 1;
    
    data = permute(data, trp_idx);

    si_trp = size(data);

    if strcmp(ftype,'bandpass')
        myFilter = design(fdesign.bandpass('N,F3dB1,F3dB2',order,fc(1),fc(2),fs),'butter');
        y=filtfilt(myFilter.sosMatrix,myFilter.ScaleValues,double(data));
    elseif strcmp(ftype,'low')
        fn=fs/2;
        wn=fc/fn;
        [B,A]=butter(order,wn,ftype);
        na=length(A);
        nb=length(B);
        nfilt = max(nb,na);   
        nfact = max(1,3*(nfilt-1));  % length of edge transients
        if  size(data,1) <= nfact
            y=double(data);
        else
            y=filtfilt(B,A,double(data));
        end
    end
    y=permute(y,trp_idx);
end
