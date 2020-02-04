%try to load the dataset
source "./help/geometry_helpers_2d.m"
source "./help/utility.m"
source "./graph_slam.m"

file = "./slam2d_range_only_initial_guess.g2o";
file_ground_truth = "./slam2d_range_only_ground_truth.g2o";

[id2state_landmark,id2state_pose,state2id_landmark,state2id_pose,associations_p_l,associations_p_p,Z_range,XR_guess,XL_guess,Zij] = create_initial_guesses(file);

[XR_ground_truth,XL_ground_truth] = load_ground_truth(file_ground_truth);

num_poses = length(state2id_pose);
num_landmarks = length(state2id_landmark);
num_iterations = 50;
damping = 0.01;
kernel_threshold = 1;


[XR, XL, chi_stats,chi_stats_p_l,chi_stats_p_p, num_inliers]=doMultiICP(XR_guess, XL_guess, Z_range, Zij,
							associations_p_l,
							associations_p_p, 
							num_poses, 
							num_landmarks, 
							num_iterations, 
							damping, 
							kernel_threshold,1);

