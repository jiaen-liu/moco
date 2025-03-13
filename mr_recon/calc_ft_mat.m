function [A,k_out,mato]=calc_ft_mat(rsamp,nk,flat,ramp,nr,apodization)
    
    if nargin==1
        para=rsamp;
        nr=para.nr;
        [rsamp, nk, flat, ramp]=cal_ramp_samp_para(para);
        apodization=0.3;
    elseif nargin==2
        apodization=nk;
        para=rsamp;
        nr=para.nr;
        [rsamp, nk, flat, ramp]=cal_ramp_samp_para(para);
    elseif nargin==5
        apodization=0.3;
    elseif nargin<5 && nargin>2
        error('*** The number of input is not supported! ***');
    end
    nr_os = nr*2;
    A=zeros(nk,nr_os);
    pos=[0:nr_os-1].'-floor(nr_os/2);
    k=[0:nk-1].';
    if flat < nk 
        % test if there are any ramp sampling
        % taken from Peter's code
        hflat= flat*0.5; %  half of flat part
        rpart= rsamp*ramp; % used part of ramp

        d= mod((hflat + rpart +0.5),1.0); % offset of samples with respect to center
        nf1= floor(hflat-d)+1; % left side of flat
        nf2= floor(flat)-nf1; % right side of flat
        nr1= floor(hflat+ rpart+ 0.5)-nf1; % left ramp
        nr2= nk- nr1- nf1- nf2; % right ramp

        for i= -nf1+1:nf2 
            k(nr1+nf1+i-1+1)= -d +i; % fill flat part
        end
        dr= nf1+ d- hflat;
        opp= -hflat- dr+ dr^2/ramp/2;
        for i= 1:nr1
            k(i)=opp-(nr1-i)+(nr1-i)^2/ramp/2;
        end
        dr=(nf2+1)-d-hflat;
        opp=hflat+dr-dr^2/ramp/2;
        for i=1:nr2
            k(nr1+nf1+nf2+i)=opp+i-1-(i-1)^2/ramp/2;
        end
        opp=flat+rpart*2-ramp*rsamp^2;
        k=k/opp*2*pi; % normalize to 2Pi
        
    else 
        k = (([0:nk-1].'-floor(nk/2.0)))/nk*2*pi; % linear sampling, no ramps involved
    end
    k_out = k;
    for i = 1:nr_os 
        A(:,i) = exp(1i*k*pos(i));
    end
    mato=A;
    A=inv(A'*A)*A';
    if apodization ~= -1
        % apodization
        k_min = min(k(:));
        k_max = max(k(:));
        w = tukeyWinNUGrid(nk, apodization, (k(:)-k_min)/(k_max-k_min));
        A=A*diag(w);
    end
end
