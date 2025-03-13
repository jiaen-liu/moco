function a=rotm2ang(m)
% rotation order is x (first), y and z
% this function is only valid for angles in the (-pi/2, pi/2) range
    si=size(m);
    ndim=numel(si);
    if ndim==2
        n=1;
    else
        n=si(3);
    end
    a=zeros(3,n);
    ax=0;
    ay=0;
    az=0;
    for i=1:n
        if m(1,1,i)^2+m(2,1,i)^2 == 0
            % cos(thetay)=0
            error('This function only works for angles in the (-pi/2, pi/2) range!');
        else
            ay=asin(-m(3,1,i));
            az=atan2(m(2,1,i),m(1,1,i));
            ax=atan2(m(3,2,i),m(3,3,i));
            a(:,i)=[ax,ay,az]';
        end
    end
end
