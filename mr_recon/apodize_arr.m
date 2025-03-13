function d=apodize_arr(d, r, dim)
  si = size(d);
  w = tukeywin(si(dim), r);
  if dim>1
      w=reshape(w,[ones(1,dim-1),si(dim)]);
  else
      w=w(:);
  end
  d=d.*w;
end
