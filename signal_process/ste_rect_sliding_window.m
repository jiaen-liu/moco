function [y1,y2]=ste_rect_sliding_window(d,sigma)
    [nr,nt]=size(d);
    dave=mean(d,1);
    ddif=d-dave;
    ddif_filter=gauss_fil_fft(ddif,sigma,1,2);
    y1=d-ddif_filter;
    y2=d-ddif;
end