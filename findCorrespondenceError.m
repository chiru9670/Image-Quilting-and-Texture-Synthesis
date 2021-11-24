%% Helper function to return correspondence error b/w target and current patch
function corresp_err = findCorrespondenceError(sel_patch,target_p,corr_type)
	if strcmp(corr_type,'luminance')
		sel_patch = rgb2hsv(sel_patch);
		target_p = rgb2hsv(target_p);
		sel_patch = sel_patch(:,:,3);
		target_p = target_p(:,:,3);
		corresp_err = rmsError(sel_patch,target_p);
    elseif strcmp(corr_type,'intensity')
		sel_patch = rgb2gray(sel_patch);
		target_p = rgb2gray(target_p);
		corresp_err = rmsError(sel_patch,target_p);
	end
end