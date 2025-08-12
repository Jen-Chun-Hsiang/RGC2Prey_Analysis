figure;
colormap(gray);
maxv = max(vol_loc0(:));
minv = min(vol_loc0(:));
for i = 1:size(vol_loc0, 1)
    imagesc(squeeze(vol_loc0(i, :, :)), [minv, maxv]);
    pause(0.5);
end
%%
figure; 
subplot(1, 3, 1);
colormap(gray)
a = squeeze(-std(vol_loc0, [], 1))';
% imagesc(a);
imagesc(a, [min(a(:)) max(a(:))+0.3*range(a(:))]);
axis off

subplot(1, 3, 2);
colormap(gray)
a = squeeze(multi_opt_sf(:, :, 250))';
a_max = max(a(:));
a_mid = median(a(:));
a_half = a_max-a_mid;
imagesc(a, [a_mid-a_half, a_max]); 
axis off
subplot(1, 3, 3);
a = (a-min(a(:)))/range(a(:));
imshow(a)

%%
figure; 
subplot(1, 3, 1);
colormap(gray)
a = squeeze(-std(vol_loc1, [], 1))';
% imagesc(a);
imagesc(a, [min(a(:)) max(a(:))+0.3*range(a(:))]);
axis off

subplot(1, 3, 2);
colormap(gray)
a = squeeze(multi_opt_sf(:, :, 120))';
a_max = max(a(:));
a_mid = median(a(:));
a_half = a_max-a_mid;
imagesc(a, [a_mid-a_half, a_max]); 
axis off
subplot(1, 3, 3);
a = (a-min(a(:)))/range(a(:));
imshow(a)