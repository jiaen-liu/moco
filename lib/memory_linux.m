function freemem=memory_linux()
    [r,w] = unix('free | grep Mem');
    stats = str2double(regexp(w, '[0-9]*', 'match'));
    memsize = stats(1)/1e6;
    % freemem = (stats(3)+stats(5))/1e6;
    freemem = stats(end)/1024^2;
end