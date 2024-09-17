function [motpar,R1,t1]=afmat2motpar(Af,source,dims,Ac,tc,varargin)
% convert affine matrix in any coordinate system to the right-hand,
% image center-orign-based coordinate
% Ac and tc is the coordinate coversion parameters from the
% to the specific one.
    p=inputParser;
    scale=[1,1,1];
    addParameter(p,'scale',scale,@isnumeric);
    p.parse(varargin{:});
    scale=p.Results.scale;
    if nargin < 4 || isempty(Ac) || isempty(tc)
        switch source
          case 'amri'
            Ac=[1,0,0;0,1,0;0,0,1];
            tc=[(dims(1)+1)/2-floor((dims(1)+1)/2);...
                (dims(2)+1)/2-floor((dims(2)+1)/2);...
                (dims(3)+1)/2-floor((dims(3)+1)/2)];
            R2=Af(1:3,1:3);
            t2=Af(1:3,4);
          case 'afni'
            Ac=[-1,0,0;0,-1,0;0,0,1];
            tc=[-(dims(1)-1)/2,-(dims(2)-1)/2,(dims(3)-1)/2].';
            R2=Af(1:3,1:3);
            t2=Af(1:3,4);
          case 'mcf'
            % this was verified
            % resds is the actual resolution
% $$$             [im_mcf,tra_mcf,~]=reg3dv(im_test,im_ds,'method','mcf','res',resds,'verbose',1);
% A uniform grid
% $$$             mot_mcf=afmat2motpar(tra_mcf(:,:,2),'mcf',[nx,ny,nz],[],[],'scale',res);
            Af=inv(Af);
            Ac=[-1,0,0;0,1,0;0,0,1];
            tc=(dims(:)-1)/2;
            R2=Af(1:3,1:3);
            t2=Af(1:3,4);
          case 'matlab'
            Af=inv(Af.');
            Ac=[0,1,0;1,0,0;0,0,1];
            tc=[(dims(2)+1)/2,(dims(1)+1)/2,(dims(3)+1)/2].';
            R2=Af(1:3,1:3);
            t2=Af(1:3,4);
          case 'rev_matlab'
            % convert into matlab's coordinate system
            Ac=[0,1,0;1,0,0;0,0,1];
            tc=-[(dims(1)+1)/2,(dims(2)+1)/2,(dims(3)+1)/2].';
            R2=Af(1:3,1:3);
            t2=Af(1:3,4);
          case 'rev_afni'
            % convert into afni's coordinate system
            Ac=[-1,0,0;0,-1,0;0,0,1];
            tc=[(dims(1)-1)/2,(dims(2)-1)/2,-(dims(3)-1)/2].';
            R2=Af(1:3,1:3);
            t2=Af(1:3,4);
        end
    elseif nargin == 4
        error('Please provide the coordinate conversion information Ac and tc!');
    end
    for i=1:3
        Ac(i,:)=Ac(i,:)*scale(i);
        tc(i)=tc(i)*scale(i);
    end
    Ac_inv=inv(Ac);
    R1=Ac_inv*R2*Ac;
    t1=Ac_inv*(t2-tc+R2*tc);

    motpar=zeros(6,1);
    motpar(1:3)=rotm2ang(R1);
    motpar(4:6)=t1;
end