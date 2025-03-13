function [b0_fit,resid,error,c_b0]=eval_sphere_fitb0(b0,mask,ord,ordn,res);
    if nargin<4
        ordn=[];
    end
    if nargin<5
        res=[1,1,1];
    end
    [nx,ny,nz]=size(mask);
    n=nx*ny*nz;
    nv=numel(b0)/n;
    x=([1:nx]-(1+nx)/2)*res(1);
    y=([1:ny]-(1+ny)/2)*res(2);
    z=([1:nz]-(1+nz)/2)*res(3);
    c_b0=sphere_harm_model_3d(b0,x,y,z,ord,mask);
    % exclude fitting to certain order
    if ~isempty(ordn)
        for i=1:numel(ordn)
            c_b0((ordn(i)+1)^2-2*ordn(i):(ordn(i)+1)^2,:)=0;
        end
    end
    [x,y,z]=ndgrid(x,y,z);
    % calculate fitted b0
    b0_fit=reshape(sphere_harm_calc_3d(x,y,z,c_b0),[nx,ny,nz,nv]);
    resid=b0-b0_fit;
    resid=reshape(resid,[n,nv]);
    resid(~mask,:)=0;
    error=a2v(std(resid(mask,:),0,1));
    resid=reshape(resid,[nx,ny,nz,nv]);
end