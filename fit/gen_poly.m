function A=gen_poly(x,y,z,ord)
    if nargin==4
        nd=3;
    elseif nargin==3
        nd=2;
        ord=z;
    elseif nargin==2
        nd=1;
        ord=y;
    end
    n=numel(x);
    poly_mat=polyns(nd+1,ord);
    size_poly_mat=size(poly_mat);
    n_poly=size_poly_mat(2);

    A=zeros(n,n_poly);
    A(:,1)=1;

    if n_poly>1
        switch nd
          case 1
            for i=2:n_poly
                A(:,i)=x.^poly_mat(2,i);
            end
          case 2
            for i=2:n_poly
                A(:,i)=x.^poly_mat(2,i).*y.^poly_mat(3,i);
            end
          case 3
            for i=2:n_poly
                A(:,i)=x.^poly_mat(2,i).*y.^poly_mat(3,i).*z.^poly_mat(4,i);
            end
        end
    end
end