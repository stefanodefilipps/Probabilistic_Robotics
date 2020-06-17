%try to load the dataset
source "./help/geometry_helpers_2d.m"
source "./help/utility.m"
source "./graph_slam.m"

file = "./slam2d_range_only_initial_guess.g2o";
file_ground_truth = "./slam2d_range_only_ground_truth.g2o";

[id2state_landmark,id2state_pose,state2id_landmark,state2id_pose,associations_p_l,associations_p_p,Z_range,XR_guess,XL_guess,Zij] = create_initial_guesses(file);

[XR_ground_truth,XL_ground_truth,Zij_ground_truth,id2state_landmark_ground_truth] = load_ground_truth(file_ground_truth);

XL_ground_truth_temp = XL_ground_truth(:,1:(size(XL_ground_truth,2)-1));
XL_ground_truth_temp(:,39:size(XL_ground_truth_temp,2)) = XL_ground_truth(:,40:end);
XR_ground_truth_odom(:,:,1) = XR_ground_truth(:,:,1);
for(i=1:size(Zij,3))
	XR_ground_truth_odom(:,:,i+1) = XR_ground_truth_odom(:,:,i)*Zij_ground_truth(:,:,i); 
end

% This is needed when i want to try to use initial guess the ground truth because I need to reorder the XL vector based on the associations
XL_guess_temp = zeros(2,size(XL_guess,2));
for(i=1:size(XL_guess,2))
	id_land = state2id_landmark(i);
	state_land = id2state_landmark_ground_truth(id_land);
	init_guess = XL_ground_truth(:,state_land);
	XL_guess_temp(:,i) = init_guess;
end

num_poses = length(state2id_pose);
num_landmarks = length(state2id_landmark);
num_iterations = 20;
damping = 0.01;
kernel_threshold = 1;


[XR, XL, chi_stats,chi_stats_p_l,chi_stats_p_p, num_inliers]=doMultiICP(XR_ground_truth, XL_guess_temp, Z_range, Zij_ground_truth,
							associations_p_l,
							associations_p_p, 
							num_poses, 
							num_landmarks, 
							num_iterations, 
							damping, 
							kernel_threshold,0);

%{
[XR, XL, chi_stats,chi_stats_p_l,chi_stats_p_p, num_inliers]=doMultiICP(XR_guess, XL_guess, Z_range, Zij,
							associations_p_l,
							associations_p_p, 
							num_poses, 
							num_landmarks, 
							num_iterations, 
							damping, 
							kernel_threshold,0);
%}

							

