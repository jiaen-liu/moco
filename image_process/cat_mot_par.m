function y=cat_mot_par(mot_par,mot_par0)
    si=size(mot_par);
    si_tmp=[6,prod(si)/6];
    mot_par=reshape(mot_par,si_tmp);
    rot_mat0=rot3d(squeeze(mot_par0(1:3)));
    y=zeros(si_tmp);
    for i=1:si_tmp(2)
        rot_mat1=rot3d(mot_par(1:3,i));
        rot_mat=rot_mat1*rot_mat0;
        tra=rot_mat1*a2v(mot_par0(4:6))+a2v(mot_par(4:6,i));
        y(4:6,i)=tra;
        y(1:3,i)=rotm2ang(rot_mat);
    end
    y=reshape(y,si);
end