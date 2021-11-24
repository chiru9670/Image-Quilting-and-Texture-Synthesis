%% MainScript
clear variables
% close all

%% Setting the color scale
my_num_of_colors = 256;
col_scale =  [0:1/(my_num_of_colors-1):1]';
my_color_scale = [col_scale,col_scale,col_scale];

%% Set to_save to 1, if you want to save the generated pictures
to_save = 0;

%% Loading the pictures
% tic;
% For GIF pictures, need to convert from index to rgb
% texture_pic = imread('data/transfer/neworange.jpg');
% texture_pic = double(texture_pic)/255.0;
[texture_pic,map] = imread('data/paper/S27.gif');
texture_pic = ind2rgb(texture_pic,map);
texture_pic = double(texture_pic);

target_pic =  imread('data/transfer/girl.jpg');
% target_pic = imresize(target_pic,0.25,'bicubic');
target_pic = double(target_pic)/255.0;
[h,w,num_chan] = size(texture_pic);
[h1,w1,num_chan] = size(target_pic);
file_name = "Result";
mid_name = "Target";
% title_name = ['Modified Pic P: ', 

%% Defining the parameters of our algorithm
patch_size = 18;
overlap_size = patch_size/3;
net_patch_size = patch_size-overlap_size;
error_tolerance = 0.1;
alph = 0.02;
num_iter = 2;
corr_type = 'intensity';
title_name = ['Modified Pic; (P: ', num2str(patch_size), ', err-tol: ', num2str(error_tolerance), ', alpha: ',num2str(alph),')']; 
%% Calculating the new generated image size 
% hnew = net_patch_size*floor(h/net_patch_size) + overlap_size;
% wnew = net_patch_size*floor(w/net_patch_size) + overlap_size;
% 
% target_pic = imresize(target_pic,[hnew wnew], 'bicubic');
% 
% modified_pic = zeros([hnew,wnew,num_chan]);
% 
% i_limit = (hnew-overlap_size)/net_patch_size;
% j_limit = (wnew-overlap_size)/net_patch_size;

tic;
%% Calculating the new generated image size 
hnew = net_patch_size*floor(h1/net_patch_size) + overlap_size;
wnew = net_patch_size*floor(w1/net_patch_size) + overlap_size;
% hnew = net_patch_size*floor(h/net_patch_size) + overlap_size;
% wnew = net_patch_size*floor(w/net_patch_size) + overlap_size;
target_pic = imresize(target_pic,[hnew wnew], 'bicubic');
target_old_pic = target_pic;

modified_pic = zeros([hnew,wnew,num_chan]);

f = waitbar(0,"Transferring Texture");
stepnum = 0;
i_limit = (hnew-overlap_size)/net_patch_size;
j_limit = (wnew-overlap_size)/net_patch_size;

%% Compute required fft
if(strcmp(corr_type, 'intensity'))
    op_gr = rgb2gray(texture_pic);
else
    op_gr = rgb2hsv(texture_pic);
    op_gr = op_gr(:,:,3);
end
    
Io_gr_fft = fft2(padarray(op_gr,[h-1 w-1],'post'));
Io_fft = fft2(padarray(texture_pic,[h-1 w-1],'post'));
Io_2 = sum(texture_pic .^ 2, 3);
op_gr_2 = op_gr.^2;
Q1 = ones(patch_size);
Q1_ext = padarray(Q1, [h-patch_size, w-patch_size], 'post');
Q1_ext_corr = xcorr2(Q1_ext, op_gr_2);
Q1_ext_corr=Q1_ext_corr(end:-1:1,end:-1:1);
Q1_ext_corr = Q1_ext_corr(h:end,w:end);

for i = 1:i_limit
	for j = 1:j_limit

		if i==1 && j==1
			target_patch = target_pic(1:patch_size,1:patch_size,:);
			modified_pic(1:patch_size,1:patch_size,:) = getFirstTransferPatch(texture_pic,target_patch,patch_size,0.0,corr_type);

		elseif i==1
            Q = ones(patch_size, overlap_size);
            Q_ext = padarray(Q, [h-patch_size, w-overlap_size], 'post');
            Q_ext_corr = xcorr2(Q_ext, Io_2);
            Q_ext_corr=Q_ext_corr(end:-1:1,end:-1:1);
            Q_ext_corr = Q_ext_corr(h:end,w:end);
			
            start_ind = net_patch_size + (j-2)*net_patch_size;
			prev_patch = modified_pic(1:patch_size,start_ind - net_patch_size + 1:start_ind - net_patch_size + patch_size,:);
			
			ref_patches = cell(1,3);
			ref_patches{1} = prev_patch;
	
			target_patch = target_pic(1:patch_size,start_ind+1:start_ind+patch_size,:);
	
			selected_patch = findClosestTransferPatch(ref_patches, target_patch, texture_pic, error_tolerance, 'vertical', overlap_size, patch_size, alph, corr_type, Io_fft, Q_ext_corr, Q1_ext_corr, Io_gr_fft );
			final_patch = minErrorBoundaryCut(ref_patches,selected_patch,overlap_size,'vertical',patch_size);
			
			modified_pic(1:patch_size,start_ind+1:start_ind+patch_size,:) = final_patch;

		elseif j==1
			start_ind = net_patch_size + (i-2)*net_patch_size;
			prev_patch = modified_pic(start_ind - net_patch_size + 1:start_ind - net_patch_size + patch_size,1:patch_size,:);
			
			ref_patches = cell(1,3);
			ref_patches{2} = prev_patch;
			
			target_patch = target_pic(start_ind+1:start_ind+patch_size,1:patch_size,:);

			selected_patch = findClosestTransferPatch(ref_patches, target_patch, texture_pic, error_tolerance, 'horizontal', overlap_size, patch_size, alph, corr_type, Io_fft, Q_ext_corr, Q1_ext_corr, Io_gr_fft );
			final_patch = minErrorBoundaryCut(ref_patches,selected_patch,overlap_size,'horizontal',patch_size);
			
			modified_pic(start_ind+1:start_ind+patch_size,1:patch_size,:) = final_patch;

		else
			left_ind = net_patch_size + (j-2)*net_patch_size;
			top_ind = net_patch_size + (i-2)*net_patch_size;
			
			left_patch = modified_pic(top_ind + 1 : top_ind + patch_size,left_ind - net_patch_size + 1:left_ind - net_patch_size + patch_size,:);
			top_patch = modified_pic(top_ind - net_patch_size + 1:top_ind - net_patch_size + patch_size,left_ind + 1:left_ind + patch_size,:);
			corner_patch = modified_pic(top_ind - net_patch_size + 1:top_ind - net_patch_size + patch_size,left_ind - net_patch_size + 1:left_ind - net_patch_size + patch_size,:);
			
			ref_patches = cell(1,3);
			ref_patches{1} = left_patch;
			ref_patches{2} = top_patch;
			ref_patches{3} = corner_patch;

			target_patch = target_pic(top_ind+1:top_ind+patch_size,left_ind+1:left_ind+patch_size,:);

			selected_patch = findClosestTransferPatch(ref_patches, target_patch, texture_pic, error_tolerance, 'both', overlap_size, patch_size, alph, corr_type, Io_fft, Q_ext_corr, Q1_ext_corr, Io_gr_fft );
			final_patch = minErrorBoundaryCut(ref_patches,selected_patch,overlap_size,'both',patch_size);

			modified_pic(top_ind+1:top_ind+patch_size,left_ind+1:left_ind+patch_size,:) = final_patch;

		end
		stepnum = stepnum + 1;
		waitbar(stepnum/(i_limit*j_limit),f,"Transferring Texture");
	end
end
 
close(f);
imwrite(modified_pic,'results/transfer/Intensity/Result.jpg');

saveFigure3(my_color_scale,texture_pic,target_old_pic,modified_pic,mid_name,title_name,file_name,1,to_save);
toc;