function print_countdown(n,i,str)
    if i==1
        fprintf(str);
    end

    if i==n
        fprintf('%d.\n', n-i);        
    else
        fprintf('%d, ', n-i);
    end
end