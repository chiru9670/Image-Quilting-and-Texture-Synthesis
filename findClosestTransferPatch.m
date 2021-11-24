%% Helper function to return a random patch of the closest patch to the given patch in the original image
function selected_patch = findClosestTransferPatch(ref_patches,target_patch,texture_pic,error_tolerance,overlap_type,overlap_size,patch_size,alph,corr_type, Io_fft, Q_ext_corr, Q1_ext_corr, Io_gr_fft);
	[h,w,num_chan] = size(texture_pic);
	num_rows = h-patch_size+1;
	num_cols = w-patch_size+1;
	error_patch = zeros([num_rows,num_cols]);
	min_error = 10000000.0;
% 	for i=1:h-patch_size+1
% 		for j=1:w-patch_size+1
% 			curr_patch = texture_pic(i:i+patch_size-1,j:j+patch_size-1,:);
% 			overlap_error = findError(curr_patch,ref_patches,overlap_type,overlap_size,patch_size);
% 			correspondence_error = 	findCorrespondenceError(curr_patch,target_patch,corr_type);
% 			total_error = alph*overlap_error + (1-alph)*correspondence_error;
% 			if total_error == 0
% 				error_patch(i,j) = 0;
% 			else
% 				error_patch(i,j) = total_error;
% 				if total_error < min_error
% 					min_error = total_error;
% 				end
% 			end			
% 		end
% 	end

    if(strcmp(overlap_type, 'vertical'))
        l_patch = ref_patches{1};
        l_patch = l_patch(:,patch_size-overlap_size+1:patch_size,:); 
        P_old = padarray(l_patch, [0 patch_size-overlap_size], 'post');
        P_norm = sum(l_patch.^2, 'all');
        P_ext = padarray(P_old, [h-patch_size, w-patch_size], 'post');
        P_ext_inv = P_ext(end:-1:1, end:-1:1,:);
        P_ext_inv_pad = padarray(P_ext_inv, [h-1, w-1], 'post');
        P_fft = fft2(P_ext_inv_pad);
        PI_ifft = ifft2(P_fft.*Io_fft);
         D = -2*PI_ifft;
        D = P_norm+ sum(D,3);
%         D = D(h:end, w:end);
%         D = D + Q_ext_corr;
        
        
    elseif strcmp(overlap_type, 'horizontal')
        l_patch = ref_patches{2};
        l_patch = l_patch(patch_size-overlap_size+1:patch_size,:,:);
        P_old = padarray(l_patch, [patch_size-overlap_size 0], 'post');
        P_norm = sum(l_patch.^2, 'all');
        P_ext = padarray(P_old, [h-patch_size, w-patch_size], 'post');
        P_ext_inv = P_ext(end:-1:1, end:-1:1,:);
        P_ext_inv_pad = padarray(P_ext_inv, [h-1, w-1], 'post');
        P_fft = fft2(P_ext_inv_pad);
        PI_ifft = ifft2(P_fft.*Io_fft);
         D = -2*PI_ifft;
        D = P_norm+ sum(D,3);
%         D = D(h:end, w:end);
%         D = D + Q_ext_corr;
        
        
    else
%         P_old = zeros(patch_size);
        l_patch = ref_patches{1};
		l_patch = l_patch(:,patch_size-overlap_size+1:patch_size,:);
        P_old_v = padarray(l_patch, [0 patch_size-overlap_size], 'post');
		t_patch = ref_patches{2};
		t_patch = t_patch(patch_size-overlap_size+1:patch_size,:,:);
		P_old_h = padarray(t_patch, [patch_size-overlap_size 0], 'post');
		c_patch = ref_patches{3};
		c_patch = c_patch(patch_size-overlap_size+1:patch_size,patch_size-overlap_size+1:patch_size,:);
        P_old_c = padarray(c_patch, [patch_size-overlap_size patch_size-overlap_size], 'post');
        P_old = P_old_v + P_old_h - P_old_c;
        P_norm = sum(P_old.^2, 'all');
        P_ext = padarray(P_old, [h-patch_size, w-patch_size], 'post');
        P_ext_inv = P_ext(end:-1:1, end:-1:1,:);
        P_ext_inv_pad = padarray(P_ext_inv, [h-1, w-1], 'post');
        P_fft = fft2(P_ext_inv_pad);
        PI_ifft = ifft2(P_fft.*Io_fft);
        D = -2*PI_ifft;
        D = P_norm+ sum(D,3);
%         D = D(h:end, w:end);
%         D = D + Q_ext_corr;
        
    end
    D = alph*D;
    target_patch = rgb2gray(target_patch);
    T_norm = sum(target_patch.^2, 'all');
    T_ext = padarray(target_patch, [h-patch_size, w-patch_size], 'post');
	T_ext_inv = T_ext(end:-1:1, end:-1:1);
    T_ext_inv_pad = padarray(T_ext_inv, [h-1, w-1],'post');
    T_fft = fft2(T_ext_inv_pad);
    TI_ifft = ifft2(T_fft.*Io_gr_fft);
    D = D + (1-alph)*(T_norm - 2*TI_ifft);
    D = D(h:end, w:end);
    D = D + alph*Q_ext_corr + (1-alph)*Q1_ext_corr;
    
    D=abs(D(1:num_rows, 1:num_cols));
    D(D<1)=1;
    min_error = min(D, [], 'all');
	min_error = min_error*(1+error_tolerance);

	[close_patches_i, close_patches_j] = find(D<min_error);
	x = randi(length(close_patches_i),1);

	start_i = close_patches_i(x);
	start_j = close_patches_j(x);

	selected_patch = texture_pic(start_i:start_i+patch_size-1, start_j:start_j+patch_size-1, :);
end