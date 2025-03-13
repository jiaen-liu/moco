function y=evalmaskbit(mdh,i)
    d=mdh.evalmask;
    d=d(:);
    d=typecast(d,'uint64');
    y=double(bitget(d,i));
end