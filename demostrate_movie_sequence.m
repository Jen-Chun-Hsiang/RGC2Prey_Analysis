id = 200;
maxv = max('sequence', [], 'all');
minv = min('sequence', [], 'all');
colormap('gray')
num_frame = size(sequence, 2);
figure; 
for i = 1:num_frame
    subplot(1, 2, 1)
    imagesc(squeeze(sequence(:, i, 1, :, :)));colorbar
    subplot(1, 2, 2)
    imagesc(squeeze(sequence(:, i, 2, :, :)));colorbar
    pause(0.05)
    sgtitle(sprintf('frame %d/%d', i, num_frame))
end