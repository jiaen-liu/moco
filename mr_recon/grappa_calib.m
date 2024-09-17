function [w,ker]=grappa_calib(dcal,par,varargin)
% $$$     p=inputParser;
    if isfield(par,'chi')
        chi=par.chi;
    else
        chi=0;
    end
% $$$     addParameter(p,'lambda',lambda,@isnumeric);
% $$$     p.parseinput(varargin{:});
% $$$     lambda=p.Results.lambda;
% dcal is full-fov calibration data
    getvar_struct(par,'ncalx','ncaly','ncalz',...
                      'nd','ry','rz','dz',...
                      'nkx','nky','nkz');
    %
    if nd==2
        [nxd,nyd,nch]=size(dcal);
        nzd=1;
    elseif nd==3
        [nxd,nyd,nzd,nch]=size(dcal);
    end
    
    nk=nkx*nky*nkz;
    nk_ch=nk*nch;
    ncal=ncalx*ncaly*ncalz;

    % check to make sure ncalx<nxd
    if (ncalx>nxd) ||...
               (ncaly>nyd) ||...
               (ncalz>nzd)
        error('*** The calibration region is too big! ***')
    end
    
    % dcal is always 3d in this program
    dcal=reshape(dcal,[nxd,nyd*nzd,nch]);
    
    ker=grappa_ker(ry,rz,dz,nky,nkz);
    nr=length(ker);
    w=zeros(nk_ch,nch,nr);
    % build calibration data matrix Adcal
    ikx=[0:nkx-1]-floor(nkx/2);
    icx_start=ceil((nxd+1)/2)-floor(ncalx/2)-1;
    icy_start=ceil((nyd+1)/2)-floor(ncaly/2)-1;
    icz_start=ceil((nzd+1)/2)-floor(ncalz/2)-1;
    Adcal=zeros(ncal,nkx,nky*nkz,nch,nr);
    bcal=zeros(ncal,nch,nr);

    
    for ic=1:ncal
        [icx,icy,icz]=ind2sub([ncalx,...
                            ncaly,...
                            ncalz],ic);
        
        
        icx=icx+icx_start;
        icy=icy+icy_start;
        icz=icz+icz_start;
        icyz=icy+(icz-1)*nyd;
        for ir=1:nr
            indyz=sub2ind([nyd,nzd],...
                          iwrapToN(icy+ker{ir}.ys,nyd),...
                          iwrapToN(icz+ker{ir}.zs,nzd));
            Adcal(ic,:,:,:,ir)=dcal(iwrapToN(icx+ikx,nxd),...
                                    indyz,:);

        end
        bcal(ic,:,:)=repmat(a2v(dcal(icx,icyz,:)),...
                            [1,nr]);
    end
    Adcal=reshape(Adcal,ncal,nk_ch,nr);
    % calculate the weights
% $$$     parfor ir=1:nr
% $$$         w(:,:,ir)=Adcal(:,:,ir)\bcal(:,:,ir);
% $$$     end
    % regularization
    parfor ir=1:nr
        Ah=Adcal(:,:,ir)';
        A=Adcal(:,:,ir);
        btmp=Ah*bcal(:,:,ir);
        M=Ah*A;
        t=trace(M);
        lambda=chi*t;
        w(:,:,ir)=(M+lambda*eye(size(M)))\btmp;
    end
end