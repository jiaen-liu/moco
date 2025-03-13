function s=ste_sense(data,mask,te,mot,coord_ste,coord_main,par,res)
% te in second
    if par.sense_necho==1
        data=data(:,:,:,:,1,:);
        nte=1;
    else
        nte=size(data,5);
        if nte~=length(te)
            error('*** The number of TE do not match! ***');
        end
    end
    if nargin<6
        coord_ste=[];
        coord_main=[];
    end
    
    si=size(data);
    nch=si(4);
    nx=si(1);
    ny=si(2);
    nz=si(3);
    if numel(si) > 5
        nim=si(6);
    else
        nim=1;
    end
    % align the images to the reference head position
    % tformtmp=affine3d(eye(4,4));
    fprintf('The number of sensitivity maps to be aligned is: ');
    % data=reshape(data,[si(1:3),nch,nim]);
    for i=1:nim
        % tformtmp.T=squeeze(tra_mat(:,:,i));
        data_tmp=data(:,:,:,:,:,i);
        for k=1:nte
            data_tmp(:,:,:,:,k)=move_matrix(data_tmp(:,:,:,:,k),...
                                            mot(:,i),...
                                            'inverse',...
                                            'res',res);
        end
        data(:,:,:,:,:,i)=data_tmp;
        fprintf('%d, ', nim-i);
    end
    fprintf('\n');
    if nte==1
        [s,ref]=sense_m(squeeze(data(:,:,:,:,1,:)),mask,eye(nch,nch),'cp_all',1);
    else
        ref_pha=zeros(nx,ny,nz,nim);
        for i=1:nim
            b0_tmp=b0_map(sum(data(:,:,:,:,:,i).*...
                              conj(data(:,:,:,:,1,i)),4),te);
            b0_tmp(~mask)=0;
            ref_pha(:,:,:,i)=exp(1i*b0_tmp*te(1)*2*pi);
                                               
        end
        [s,ref]=sense_m(squeeze(data(:,:,:,:,1,:)),mask,eye(nch,nch),'ref_pha_in',ref_pha);
    end
    % ----------------------------------------------- %
    % Might need to normalize sensivity phase to the refernece
    % future work
    % ----------------------------------------------- %
    % calculate B0 field from each pair of data contrasts
% $$$     b0=zeros(si(1),si(2),si(3),nim);
% $$$ 
% $$$     data_b0=squeeze(sum(data(:,:,:,:,2,:).*conj(data(:,:,:,:,1,:)),4));
% $$$     for i=1:nim
% $$$         b0(:,:,:,i)=unwrapper_3d_mask(angle(data_b0(:,:,:,i)),...
% $$$                                       ones(si(1:3)))/2/pi/(te(2)-te(1));
% $$$     end
% $$$     b0(mask)=0;
% $$$     % normalize sensitivity phase by the reference image
% $$$     s.sensit=s.sensit./...
% $$$              exp(1i*2*pi*reshape(b0,[si(1),si(2),si(3),1,nim])*te(1)).*...
% $$$              exp(1i*angle(reshape(s.ref,[si(1),si(2),si(3),1, ...
% $$$                         nim])));
    % ----------------------------------------------- %
    if isempty(coord_main)
        return;
    end
    % interpolate the sensitiivty to resolution of the main
    % acquisition
    si_main=size(coord_main);
    nx_main=si_main(1);
    ny_main=si_main(2);
    nz_main=si_main(3);
    sensit_main=zeros(nx_main,ny_main,nz_main,nch,nim);
    fprintf('The number of sensitivity maps left is: ');
    parfor i=1:nim
        sensit_main(:,:,:,:,i)=interp3_nmat(coord_ste,coord_main,...
                                            s(:,:,:,:,i));
        fprintf('%d, ', nim-i);
    end
    fprintf('\n');
    s=sensit_main;
end