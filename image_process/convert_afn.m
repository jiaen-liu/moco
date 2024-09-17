function afn_o=convert_afn(afn_i,R,A,d)
    R0=afn_i(1:3,1:3);
    T0=afn_i(1:3,4);
    if numel(A)==3
        A=diag(A);
    end
    d=d(:);
    Ainv=inv(A);
    Rinv=inv(R);
    R1=Ainv*Rinv*R0*R*A;
    T1=Ainv*Rinv*T0-Ainv*Rinv*d+Ainv*Rinv*R0*d;
    afn_o=[[R1,T1];[0,0,0,1]];
end