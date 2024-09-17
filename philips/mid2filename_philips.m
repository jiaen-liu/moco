function f=mid2filename_philips(mid)
    if ischar(mid)
        f=mid;
        return;
    end
    if ~isnumeric(mid)
        error('*** The data ID should be a number or string! ***');
    end
    %
    if ismatrix(mid) && numel(mid)==2
        str=num2str(mid(1),'_%d_');
        str=num2str(mid(2),[str,'%d_']);
        filter=['*',str,'*','.raw'];
        f=get_file_filter('.',filter);
        return;
    end
    if mod(mid,1)==0
        str=num2str(floor(mid));
        % str(str=='.')='_';
        str=['_',str,'_'];
        filter=['*',str,'*','.raw'];
        f=get_file_filter('.',filter);
        if isstring(f) && length(f)>1
            warning('*** More than one files have the same MID! Trying fractional MID again ***');
            str=num2str(floor(mid)+0.2);
            str(str=='.')='_';
            str=['_',str,'_'];
            filter=['*',str,'*','.raw'];
            f=get_file_filter('.',filter);
        end
    else
        str=num2str(mid);
        str(str=='.')='_';
        str=['_',str,'_'];
        filter=['*',str,'*','.raw'];
        f=get_file_filter('.',filter);
    end
    if isempty(f)
        error('*** The requested file does not exist! ***');
    end
end