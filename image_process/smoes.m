function y=smoes(im,frac,mask,varargin)
    fft_based=0;
    nvar=length(varargin);
    if nvar>0
        fft_based=varargin{1};
    end
    si=size(im);
    nx=si(1);
    ny=si(2);
    if numel(si)<3
        nz=1;
        nch=1;
    elseif numel(si)<4
        nz=si(3);
        nch=1;
    elseif numel(si)<5
        nz=si(3);
        nch=si(4);
    end
    y=zeros(nx,ny,nz,nch);
    [nx,ny,nz,nch]=size(im);
    nxy=nx*ny;
    sig=[nx,ny]/frac;
    sk=round_odd(sig*6+1);
    sigf=0.5/pi*frac;
    h=ndgauss(sk,sig);
    hf=ndgauss([nx,ny],[sigf,sigf]);
    hf=hf/max(abs(hf(:)));
    for iz=1:nz
        if isempty(mask)
            mb_inv=ones(nx,ny);
        else
            if ~fft_based
                mb=imfilter(mask(:,:,iz),h);
            else
                mb=abs(fftmr(fftmr(mask(:,:,iz),-1,[1,2]).*hf,1,[1,2]));
            end
            mb_inv=zeros(nx,ny);
            idx=find(mb>0.2);
            if numel(idx)>0
                mb_inv(idx)=1./mb(idx);
            end
        end
        for ich=1:nch
            if ~fft_based
                y(:,:,iz,ich)=imfilter(im(:,:,iz,ich).*mask(:,:,iz),h).*mb_inv;
            else
                y(:,:,iz,ich)=fftmr(fftmr(im(:,:,iz,ich).*mask(:,:,iz),-1,[1,2]).*hf,...
                                    1,[1,2]).*mb_inv;
            end
        end
    end
end