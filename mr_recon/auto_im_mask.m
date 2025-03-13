function mask=auto_im_mask(mag)
[~,thrd]=mask1d(mag(:),0.3,0.1);
mask=mag>thrd;
end