%% Texture transfer
clear variables
% close all

%% Load required pics
% tic;
% For GIF, convert index to rgb.

target =  imread('data/transfer/bill-big.jpg');
target = double(target)/255.0;

texture = imread('data/transfer/rice.jpg');
texture = double(texture)/255.0;
% [texture,map] = imread('data/paper/S27.gif');
% texture = ind2rgb(texture,map);
% texture = double(texture);

[h,w,num_channels] = size(texture);
[h1,w1,num_channels] = size(target);
file_name = "Result";
mid_name = "Target";


%% Defining params
err_tol = 0.1;
alpha = 0.3;

patch_dim = 18;
overlap_size = patch_dim/6;
net_patch_dim = patch_dim-overlap_size;

corr_type = 'intensity';
title_name = ['Modified Pic; (P: ', num2str(patch_dim), ', err-tol: ', num2str(err_tol), ', alpha: ',num2str(alpha),')']; 


tic;
%% size of the modified img
hnew = net_patch_dim*floor(h1/net_patch_dim) + overlap_size;
wnew = net_patch_dim*floor(w1/net_patch_dim) + overlap_size;

target = imresize(target,[hnew wnew], 'bicubic');
target_old = target;

result = zeros([hnew,wnew,num_channels]);


lim_x = (hnew-overlap_size)/net_patch_dim;
lim_y = (wnew-overlap_size)/net_patch_dim;

%% Compute required fft
if(strcmp(corr_type, 'intensity'))
    op_gr = rgb2gray(texture);
else
    op_gr = rgb2hsv(texture);
    op_gr = op_gr(:,:,3);
end
    
Io_gr_fft = fft2(padarray(op_gr,[h-1 w-1],'post'));
Io_fft = fft2(padarray(texture,[h-1 w-1],'post'));
Io_2 = sum(texture .^ 2, 3);
op_gr_2 = op_gr.^2;
Q1 = ones(patch_dim);
Q1_ext = padarray(Q1, [h-patch_dim, w-patch_dim], 'post');
Q1_ext_corr = xcorr2(Q1_ext, op_gr_2);
Q1_ext_corr=Q1_ext_corr(end:-1:1,end:-1:1);
Q1_ext_corr = Q1_ext_corr(h:end,w:end);

for i = 1:lim_x
	for j = 1:lim_y

		if i==1 && j==1
			target_patch = target(1:patch_dim,1:patch_dim,:);
			result(1:patch_dim,1:patch_dim,:) = getFirstTransferPatch(texture,target_patch,patch_dim,0.0,corr_type);

		elseif i==1
            Q = ones(patch_dim, overlap_size);
            Q_ext = padarray(Q, [h-patch_dim, w-overlap_size], 'post');
            Q_ext_corr = xcorr2(Q_ext, Io_2);
            Q_ext_corr=Q_ext_corr(end:-1:1,end:-1:1);
            Q_ext_corr = Q_ext_corr(h:end,w:end);
			
            start_ind = (j-1)*net_patch_dim;
			prev_patch = result(1:patch_dim,start_ind - net_patch_dim + 1:start_ind - net_patch_dim + patch_dim,:);
			
			ref_patches = cell(1,3);
			ref_patches{1} = prev_patch;
	
			target_patch = target(1:patch_dim,start_ind+1:start_ind+patch_dim,:);
	
			selected_patch = findClosestTransferPatch(ref_patches, target_patch, texture, err_tol, 'vertical', overlap_size, patch_dim, alpha, corr_type, Io_fft, Q_ext_corr, Q1_ext_corr, Io_gr_fft );
			final_patch = minErrorBoundaryCut(ref_patches,selected_patch,overlap_size,'vertical',patch_dim);
			
			result(1:patch_dim,start_ind+1:start_ind+patch_dim,:) = final_patch;

		elseif j==1
			Q = ones(overlap_size, patch_dim);
            Q_ext = padarray(Q, [h-overlap_size w-patch_dim], 'post');
            Q_ext_corr = xcorr2(Q_ext, Io_2);
            Q_ext_corr=Q_ext_corr(end:-1:1,end:-1:1);
            Q_ext_corr = Q_ext_corr(h:end, w:end);
            
            start_ind = (i-1)*net_patch_dim;
			prev_patch = result(start_ind - net_patch_dim + 1:start_ind - net_patch_dim + patch_dim,1:patch_dim,:);
			
			ref_patches = cell(1,3);
			ref_patches{2} = prev_patch;
			
			target_patch = target(start_ind+1:start_ind+patch_dim,1:patch_dim,:);

			selected_patch = findClosestTransferPatch(ref_patches, target_patch, texture, err_tol, 'horizontal', overlap_size, patch_dim, alpha, corr_type, Io_fft, Q_ext_corr, Q1_ext_corr, Io_gr_fft );
			final_patch = minErrorBoundaryCut(ref_patches,selected_patch,overlap_size,'horizontal',patch_dim);
			
			result(start_ind+1:start_ind+patch_dim,1:patch_dim,:) = final_patch;

        else
            Q_ext = zeros(h,w);
            left_ind = (j-1)*net_patch_dim;
			top_ind = (i-1)*net_patch_dim;
			
			left_patch = result(top_ind + 1 : top_ind + patch_dim,left_ind - net_patch_dim + 1:left_ind - net_patch_dim + patch_dim,:);
			top_patch = result(top_ind - net_patch_dim + 1:top_ind - net_patch_dim + patch_dim,left_ind + 1:left_ind + patch_dim,:);
			corner_patch = result(top_ind - net_patch_dim + 1:top_ind - net_patch_dim + patch_dim,left_ind - net_patch_dim + 1:left_ind - net_patch_dim + patch_dim,:);
			
            Q_ext(1:patch_dim, 1:overlap_size) = ones(patch_dim, overlap_size);
            Q_ext(1:overlap_size, 1:patch_dim) = ones(overlap_size,patch_dim);
            Q_ext_corr = xcorr2(Q_ext, Io_2);
            Q_ext_corr=Q_ext_corr(end:-1:1,end:-1:1);
            Q_ext_corr = Q_ext_corr(h:end, w:end);
            
			ref_patches = cell(1,3);
			ref_patches{1} = left_patch;
			ref_patches{2} = top_patch;
			ref_patches{3} = corner_patch;

			target_patch = target(top_ind+1:top_ind+patch_dim,left_ind+1:left_ind+patch_dim,:);

			selected_patch = findClosestTransferPatch(ref_patches, target_patch, texture, err_tol, 'both', overlap_size, patch_dim, alpha, corr_type, Io_fft, Q_ext_corr, Q1_ext_corr, Io_gr_fft );
			final_patch = minErrorBoundaryCut(ref_patches,selected_patch,overlap_size,'both',patch_dim);

			result(top_ind+1:top_ind+patch_dim,left_ind+1:left_ind+patch_dim,:) = final_patch;

		end
% 		stepnum = stepnum + 1;
% 		waitbar(stepnum/(i_limit*j_limit),f,"Transferring Texture");
	end
end


fig = figure; 
% colormap(my_color_scale);

colormap jet

subplot(1,3,1), imagesc(texture), title('Original Image'), colorbar, daspect([1 1 1]), axis tight;
subplot(1,3,2), imagesc(target_old), title(mid_name), colorbar, daspect([1 1 1]), axis tight;
subplot(1,3,3), imagesc(result), title(title_name), colorbar, daspect([1 1 1]), axis tight;


toc;

%%
% fn = 'C:\Users\rbhol\Desktop\Image-Quilting-and-Texture-Synthesis-master';
%  figlist = findobj(allchild(0),'flat','Type','figure');
%  for ifig=1:length(figlist)
%      fh=figlist(ifig);
%      fna=num2str(get(fh,'Number'));
%      set(0,'CurrentFigure',fh);
%      saveas(fh,fullfile(fn,[fna,'.jpeg']));
%  end
