function [scores,coeffs,latent,tsq,b0_fit] = ste_pca_b0(b0,mask,n_PCA_Modes,t_echo)
%Function to perform PCA and return its results
%scores are the PCA modes
%coeffs are the weights of said modes
%latent: are the principal commponent variances
%tsq: a measure of the variation explained by a given component, can be
%used in the future to automatically determine the number of components
%b0_fit: the fited b0
[nx,ny,nz]=size(mask);
n=nx*ny*nz;
nv=numel(b0)/n;

%Perform PCA
edB = exp(1i*2*pi*(b0.*mask)*t_echo*1e-3);
edBMat = reshape(edB,[n,nv]);
[coeff,score,latent,tsq] = pca(edBMat,'Centered',false);
coeffs = coeff(:,1:(n_PCA_Modes));
scores = score(:,1:(n_PCA_Modes));


%get b0_fit
edB_fit = scores*coeffs';
edB_fit = reshape(edB_fit,[nx,ny,nz,nv]);
b0_fit = angle(edB_fit)/(2*pi*t_echo*1e-3);

b0_fit = reshape(b0_fit,[nx,ny,nz,nv]).*mask;

%Reshape scores
scores = reshape(scores,[nx,ny,nz,n_PCA_Modes]);


