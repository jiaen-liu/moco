function A=covSiem(mid)
% get noise data
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
    % Matlab's cov is the conjugate of the 
    % standard definition
    A=conj(cov(n));
end