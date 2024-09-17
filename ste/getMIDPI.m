function y=getMIDPI(mid)
    if ~isstruct(mid)
        fnData=['mid' int2str(mid) '.steref4recon.svd'];
        data=read_data(rp(fnData));
    else
        data=mid;
    end
    if ~isfield(data,'mid_pimg')
        y=0;
        return;
    end
    nPI=length(data.mid_pimg);
    for i=1:nPI
        p=getfield(data.para_pimg,['mid',int2str(data.mid_pimg(i))]);
        if p.dimen==2
            y=double(data.mid_pimg(i));
            return;
        end
    end
    if nPI==1 && field_true(data,'sense_philips')
        y=data.mid_pimg;
    end
% $$$     for i=1:nPI
% $$$         p=getfield(data.para_pimg,['mid',int2str(data.mid_pimg(i))]);
% $$$         if p.dimen==3
% $$$             y=double(data.mid_pimg(i));
% $$$             return;
% $$$         end
% $$$     end
    
end
