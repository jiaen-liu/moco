function f=mid2filename(mid)
    filter=['meas_MID*',num2str(mid),'_*','.dat'];
    ntrial=5;
    name_head='meas_MID';
    for i=1:ntrial
        filter=[name_head,num2str(mid),'_*','.dat'];
        f=get_file_filter('.',filter);
        if isstring(f) && length(f)>1
            error('*** More than one files have the same MID! ***');
        end
        if ~isempty(f)
            break;
        end
        name_head=[name_head,'0'];
    end
end