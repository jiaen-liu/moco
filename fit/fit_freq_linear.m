function f=fit_freq_linear(data,t)
    nt=length(data);
    if isreal(data)
        data=exp(1i*data);
    end
    if numel(t)==1
        t=[0:nt-1]'*t;
    end
    t=col(t);
    % fit the difference of the middle echoes
    nh=floor(nt/2);
    f1=angle(data(nh+1)/data(nh))/2/pi/(t(nh+1)-t(nh));
    % fit the other echoes
    phase=angle(data/data(1));
    t=t-t(1);
    dp_fit1=t*f1*2.0*pi;
    dp_dif=wrapToPi(phase-dp_fit1);
    kd=[t,ones(nt,1)]\dp_dif;
    f=f1+kd(1)/2/pi;
end