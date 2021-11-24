%% Helper function to return the random first patch
function selected_patch = getFirstTransferPatch(texture,target_patch,patch_dim,err_tol,corr_type)
	[h,w,num_chan] = size(texture);
	rows = h-patch_dim+1;
	cols = w-patch_dim+1;
	D = zeros([rows,cols]);
	min_err = 10000000.0;
	for i=1:h-patch_dim+1
		for j=1:w-patch_dim+1
			curr_patch = texture(i:i+patch_dim-1,j:j+patch_dim-1,:);
			corresp_err = findCorrespondenceError(curr_patch,target_patch,corr_type);
			total_err = corresp_err;
			if total_err == 0
				D(i,j) = 0;
			else
				D(i,j) = total_err;
				if total_err < min_err
					min_err = total_err;
				end
			end			
		end
	end
	
	min_err = min_err*(1+err_tol);

	[res_patches_i, res_patches_j] = find(D<=min_err);
	
    x = randi(length(res_patches_i),1);

	start_i = res_patches_i(x);
	start_j = res_patches_j(x);

	selected_patch = texture(start_i:start_i+patch_dim-1, start_j:start_j+patch_dim-1, :);
end