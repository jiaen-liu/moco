function afmat=affine_forward(mot_par)
    mot_par=squeeze(mot_par);
    nv=size(mot_par,2);
    afmat=zeros(4,4,nv);
    afmat(4,4,:)=1;
    for i=1:nv
        afmat(1:3,1:3,i)=rot3d(mot_par(1:3,i));
        afmat(1:3,4,i)=mot_par(4:6,i);
    end
end