function y=save_afnmat_fsl(folder,A)
    [n1,n2,v]=size(A);
    if n1~=4 || n2~=4
        error('*** The affine matrix should have 4 rows and 4 columns! ***');
    end
    % check if folder exist
    if exist(folder)
        % remove the folder
        rmdir(folder,'s');
    end
    % make a folder
    mkdir(folder);
    for i=1:v
        fn=fullfile(folder,['MAT_' sprintf('%05d',i-1)]);
        dlmwrite(fn,A(:,:,i),'delimiter','\t');
    end
end