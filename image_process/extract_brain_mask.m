function mask_brain=extract_brain_mask(mask)
si=size(mask);
ndim=numel(si);
mask_brain=false(si);
if ndim == 2
    ns=1;
else
    ns=si(end);
end
for i=1:ns
    
    % Extract the two largest blobs, which will either be the skull and brain,
    % or the skull/brain (if they are connected) and small noise blob.
    binaryImage = bwareafilt(mask(:,:,i), 2);		% Extract 2 largest blobs.
                                                    % Erode it a little with imdilate().
    binaryImage = imopen(binaryImage, true(5));
    % Now brain should be disconnected from skull, if it ever was.
    % So extract the brain only - it's the largest blob.
    binaryImage = bwareafilt(binaryImage, 1);		% Extract largest blob.
                                                    % Fill any holes in the brain.
    binaryImage = imfill(binaryImage, 'holes');
    % Dilate mask out a bit in case we've chopped out a little bit of brain.
    % binaryImage = imdilate(binaryImage, true(3));
    mask_brain(:,:,i)=binaryImage;
end
end