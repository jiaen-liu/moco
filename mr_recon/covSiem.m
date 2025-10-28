function A=covSiem(mid,unity_diag)
% get noise data
    if nargin<2
        unity_diag=0;
    end
    if all(size(mid)>1)
        n=mid;
        [~,nch]=size(n);
    else
        fn=['mid',num2str(mid),'_noise_raw.mat'];
        if exist(fn)
            load(fn);
        else
            n=readSortSiem(mid,'noise',1);
            save(fn,'n');
        end
        % trim the head and tail
        frac=0.9;
        [nx,necho,ns,nch]=size(n);
        nedge=floor(nx*(1-frac)/2);
        nxn=nx-2*nedge;
        n=n(nedge+1:end-nedge,:,:,:);
        n=reshape(n,[nxn*necho*ns,nch]);
    end
    % Matlab's cov is the conjugate of the 
    % standard definition
    A=conj(cov(n));
    if unity_diag~=0
        n=n./reshape(abs(A(idx_diag(A))).^0.5,[1,nch]);
        A=conj(cov(n));
    end
end
