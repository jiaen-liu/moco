function y=read_pe_files(filename)
  if ~exist(filename)
     error('*** The file does not exist! ***');
  end
  finfo = dir(filename);
  nbytes = finfo.bytes;
  offset = 0;
  f=fopen(filename,'r');
  header_flt = zeros(6,1,'single');
  header_flt=typecast(fread(f,size_of(header_flt),'uint8=>uint8'),'single');
  offset = offset+6*4;
  header_l = zeros(2,1,'int32');
  header_l=typecast(fread(f,size_of(header_l),'uint8=>uint8'),'int32');
  offset = offset+2*4;
  nkyz = header_l(1);
  kyz = zeros(2, nkyz,'single');
  kyz=reshape(typecast(fread(f,size_of(kyz),'uint8=>uint8'),'single'),...
              [2,nkyz]);
  offset = offset+2*nkyz*4;
  if nbytes > offset
     parameters = zeros((nbytes-offset)/4,1,'single');
     parameters=typecast(fread(f,size_of(parameters),'uint8=>uint8'),'single');
  end
  fclose(f);
  if nbytes > offset
      y=struct('hf',header_flt,...
               'hl',header_l,...
               'k',kyz,...
               'r',parameters(1),...
               'sp',parameters(2),...
               'dkz',parameters(3),...
               'blipless',parameters(4),...
               'ncontr',parameters(5),...
               'y_cyc',parameters(6),...
               'z_cyc',parameters(7));
  else
      y=struct('hf',header_flt,...
           'hl',header_l,...
           'k',kyz);
  end
end
