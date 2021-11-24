%% MainScript
clear variables
% close all

%% Setting the color scale
my_num_of_colors = 256;
col_scale =  [0:1/(my_num_of_colors-1):1]';
my_color_scale = [col_scale,col_scale,col_scale];



%% For GIF pictures, convert from index to rgb
inp = 'cans_sc';    
folder = 'paper/';
op_n = strcat(inp,'.jpg');
input_file = strcat(inp,'.gif');

[texture,map] = imread(strcat('data/',folder,input_file));
texture = ind2rgb(texture,map);
original_pic = double(texture);

%% For JPEG pictures

% inp = 'green_plum';
% folder = 'own/';
% op_n = strcat(inp,'.jpg');
% input_file = strcat(inp,'.jpg');
% textureP = imread(strcat('data/',folder,input_file));
% original_pic = double(textureP)/255.0;
% 
[h,w,num_chan] = size(original_pic);


%% Defining the params
err_tol = 0.1;

patch_dim = 36;
overlap_size = patch_dim/6;
net_patch_dim = patch_dim-overlap_size;

title_name = ['Quilted Pic; (P: ', num2str(patch_dim), ', err-tol: ', num2str(err_tol), ')'];

tic;
%% Calculating the new generated image size 
h_new = 2*net_patch_dim*floor(h/net_patch_dim) + overlap_size;
w_new = 2*net_patch_dim*floor(w/net_patch_dim) + overlap_size;

Quilt_pic = zeros([h_new,w_new,num_chan]);

lim_x = (h_new-overlap_size)/net_patch_dim;
lim_y = (w_new-overlap_size)/net_patch_dim;

%% compute required fft
Io_fft = fft2(padarray(original_pic,[h-1 w-1],'post'));
Io_2 = sum(original_pic .^ 2, 3);

for i = 1:lim_x
	for j = 1:lim_y

		if i==1 && j==1
			Quilt_pic(1:patch_dim,1:patch_dim,:) = getRandomPatch(original_pic,patch_dim);

		elseif i==1
			Q = ones(patch_dim, overlap_size);
            Q_ext = padarray(Q, [h-patch_dim, w-overlap_size], 'post');
            Q_ext_corr = xcorr2(Q_ext, Io_2);
            Q_ext_corr=Q_ext_corr(end:-1:1,end:-1:1);
            Q_ext_corr = Q_ext_corr(h:end,w:end);
            
            start_ind = (j-1)*net_patch_dim;
			prev_patch = Quilt_pic(1:patch_dim,start_ind - net_patch_dim + 1:start_ind - net_patch_dim + patch_dim,:);
			
			prevs_p = cell(1,3);
			prevs_p{1} = prev_patch;
			
			selected_patch = findClosestPatch(prevs_p, original_pic, err_tol, 'vertical', overlap_size, patch_dim, Io_fft, Q_ext_corr);
			final_patch = minErrorBoundaryCut(prevs_p,selected_patch,overlap_size,'vertical',patch_dim);
			
			Quilt_pic(1:patch_dim,start_ind+1:start_ind+patch_dim,:) = final_patch;

		elseif j==1
			Q = ones(overlap_size, patch_dim);
            Q_ext = padarray(Q, [h-overlap_size w-patch_dim], 'post');
            Q_ext_corr = xcorr2(Q_ext, Io_2);
            Q_ext_corr=Q_ext_corr(end:-1:1,end:-1:1);
            Q_ext_corr = Q_ext_corr(h:end, w:end);
            
            start_ind = (i-1)*net_patch_dim;
			prev_patch = Quilt_pic(start_ind - net_patch_dim + 1:start_ind - net_patch_dim + patch_dim,1:patch_dim,:);
			
			prevs_p = cell(1,3);
			prevs_p{2} = prev_patch;
			
			selected_patch = findClosestPatch(prevs_p, original_pic, err_tol, 'horizontal', overlap_size, patch_dim, Io_fft, Q_ext_corr);
			final_patch = minErrorBoundaryCut(prevs_p,selected_patch,overlap_size,'horizontal',patch_dim);
			
			Quilt_pic(start_ind+1:start_ind+patch_dim,1:patch_dim,:) = final_patch;

        else
            Q_ext = zeros(h,w);
            left_ind = (j-1)*net_patch_dim;
			top_ind = (i-1)*net_patch_dim;
			
			left_patch = Quilt_pic(top_ind + 1 : top_ind + patch_dim,left_ind - net_patch_dim + 1:left_ind - net_patch_dim + patch_dim,:);
			top_patch = Quilt_pic(top_ind - net_patch_dim + 1:top_ind - net_patch_dim + patch_dim,left_ind + 1:left_ind + patch_dim,:);
			corner_patch = Quilt_pic(top_ind - net_patch_dim + 1:top_ind - net_patch_dim + patch_dim,left_ind - net_patch_dim + 1:left_ind - net_patch_dim + patch_dim,:);
			
            Q_ext(1:patch_dim, 1:overlap_size) = ones(patch_dim, overlap_size);
            Q_ext(1:overlap_size, 1:patch_dim) = ones(overlap_size,patch_dim);
            Q_ext_corr = xcorr2(Q_ext, Io_2);
            Q_ext_corr=Q_ext_corr(end:-1:1,end:-1:1);
            Q_ext_corr = Q_ext_corr(h:end, w:end);
			
            prevs_p = cell(1,3);
			prevs_p{1} = left_patch;
			prevs_p{2} = top_patch;
			prevs_p{3} = corner_patch;

			selected_patch = findClosestPatch(prevs_p, original_pic, err_tol, 'both', overlap_size, patch_dim, Io_fft, Q_ext_corr);
			final_patch = minErrorBoundaryCut(prevs_p,selected_patch,overlap_size,'both',patch_dim);

			Quilt_pic(top_ind+1:top_ind+patch_dim,left_ind+1:left_ind+patch_dim,:) = final_patch;

		end
			
		
	end
end
toc;

fig = figure; colormap(my_color_scale);
colormap jet;
subplot(1,2,1), imagesc(original_pic), title('Original Image'), colorbar, daspect([1 1 1]), axis tight;
subplot(1,2,2), imagesc(Quilt_pic), title(title_name), colorbar, daspect([1 1 1]), axis tight;
impixelinfo();
