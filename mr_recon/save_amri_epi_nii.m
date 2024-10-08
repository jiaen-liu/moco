function status=save_amri_epi_nii(path,dat,para,prefix,min_inplane_rot)
status=1;
if nargin<5
    min_inplane_rot=1;
end
% save nifti
if size(dat,3)>para.n_partitions_nos
    dat=dat(:,:,...
            idx_truncate(para.n_partitions,...
                         para.n_partitions_nos),...
            :,:);
end
% save mean magnitude
fn=fullfile(path,[prefix,...
                  '_mag_echo_ave.nii.gz']);
siem_to_nifti(fn,single(mean(abs(dat).^2,4).^0.5),...
              para,1,1);
fn_odd=fullfile(path,[prefix,...
                  '_mag_echo_odd.nii.gz']);
siem_to_nifti(fn_odd,single(mean(abs(dat(:,:,:,1:2:end)).^2,4).^0.5),...
              para,1,1);
if size(dat,4)>1
    fn_even=fullfile(path,[prefix,...
                        '_mag_echo_even.nii.gz']);
    siem_to_nifti(fn_even,single(mean(abs(dat(:,:,:,2:2:end)).^2,4).^0.5),...
                            para,1,1);
end
% save magnitude and phase of each echo
for k=1:size(dat,4)
    fnm=fullfile(path,[prefix,...
                       '_mag_e',num2str(k),...
                       '.nii.gz']);
    fnp=fullfile(path,[prefix,...
                       '_pha_e',num2str(k),...
                       '.nii.gz']);
    siem_to_nifti(fnm,single(abs(dat(:,:,:,k,:))),...
                  para,1,min_inplane_rot);
    siem_to_nifti(fnp,single(angle(dat(:,:,:,k,:))),...
                  para,1,min_inplane_rot);
end
status=0;
end