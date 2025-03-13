function y=read_fitrec(fn,dims)
    nd=20;
    y=read_raw(fn,[],'double');
    nv=numel(y)/nd;
    if mod(nv,1) ~= 0
        error('The data is truncated!');
    end
    y=reshape(y,[nd,nv]);
    y=y([16,15,14,1,2,3],:);
    for i=1:nv
        am=rot3d(0,0,y(3,i))*...
           rot3d(0,y(2,i),0)*...
           rot3d(y(1,i),0,0);
        am=[am -squeeze(y(4:6,i))];
        am=[am;[0,0,0,1]];
        y(:,i)=afmat2motpar(am,'amri',dims);
    end
end