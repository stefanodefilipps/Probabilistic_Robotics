source "./help/geometry_helpers_2d.m"

% This is the function that parse the file containing the initial guesses on robot poses and odometry measurement and range measurement and build 
% the correct structures in order to launch the LS solver

function [id2state_landmark,id2state_pose,state2id_landmark,state2id_pose,associations_p_l,associations_p_p,Z_range,XR_guess,XL_guess,Zij] = create_initial_guesses(file) 
	f = fopen(file);
	initial_guess_poses = [];
	% vector of association for the range measurement
	associations_p_l = [];
	% vector of association for the odometry measurement
	associations_p_p = [];
	Z_range = [];
	%Xij will contain the relative displacement between position i and j so i build a kind of virtual measurement
	Zij = zeros(3,3,1);
	flag = 0;
	%bookkeeping: to and from mapping between robot pose (x,y, theta) and state vector and landmark indices (i) and state vector and viceversa
	id2state_landmark = [];
	id2state_pose = [];
	state2id_landmark = [];
	state2id_pose = [];
	land_meas_counter = [];
	l = fgetl(f);
	while ischar(l)
		d = strsplit(l);
		header = d{1};
		%get all initial guesses for the robot poses
		if strcmp(header,"VERTEX_SE2") == 1
			state2id_pose(end+1) = str2double(d{2});
			id2state_pose(str2double(d{2})) = length(state2id_pose);
			initial_guess_poses(end+1,1) = str2double(d{3});
			initial_guess_poses(end,2) = str2double(d{4});
			initial_guess_poses(end,3) = str2double(d{5});
		end
		if strcmp(header,"EDGE_SE2") == 1
			% when we are parsing the odometry measurements then we need to build the virtual measurement between pose j and pose i
			% which will be an homogenous transformation matrix representing the relative position of pose j wrt pose i
			id_pose_i = str2double(d{2});
			id_pose_j = str2double(d{3});
			dx_r = str2double(d{4});
			dy_r = str2double(d{5});
			dtheta_r = str2double(d{6});
			dXij = v2t([dx_r;dy_r;dtheta_r]);
			state_pose_i = id2state_pose(id_pose_i);
			state_pose_j = id2state_pose(id_pose_j);
			associations_p_p(:,end+1) = [state_pose_i;state_pose_j];
			if flag == 0
				Zij(:,:,1) = dXij;
				flag = 1;
			else
				Zij(:,:,end+1) = dXij;
			end
		end
		if strcmp(header,"EDGE_RANGE_SE2_XY") == 1
			% When I am parsin the measurements I only retain the landmarks for which I have at least 3 range measurements, otherwise I can't compute
			% a correct initialization of the landamrks.
			state_pose = id2state_pose(str2double(d{2}));
			id_landmark = str2double(d{3});
			% I need to put in the states only the landmark for which I have at least 3 measurements
			if length(land_meas_counter) < id_landmark 
				% In this case I have never encountered the landmark so far so I initialize the counter to 1
				land_meas_counter(id_landmark) = 1;
			elseif land_meas_counter(id_landmark) == 0
				% In this case i have already encountered a landmark whose id was bigger than current landmark_id and so land_meas_counter(id_landmark) was
				% set to 0 and I need to initialize it to 1
				land_meas_counter(id_landmark) = 1;
			else
				% I have already encountered this landmark and I simply increase the counter
				land_meas_counter(id_landmark) = land_meas_counter(id_landmark) + 1;
			end
			%if I have at least 3 measurement for a landmark then I add the landmark in the state
			if land_meas_counter(id_landmark) == 3
				state2id_landmark(end+1) = id_landmark;
				id2state_landmark(id_landmark) = length(state2id_landmark);
			end
			range_meas = str2double(d{4});
			associations_p_l(:,end+1) = [state_pose;id_landmark];
			Z_range(end+1) = range_meas;
		end
		l = fgetl(f);
	end
	temp_ass = [];
	temp_meas = [];

	% from the associations and measurement I need only to retrieve the measurements that are associated to landmark in the actual state

	for (i=1:columns(associations_p_l))
		id_landmark = associations_p_l(2,i);
		if any(state2id_landmark == id_landmark) == 1
			state_pose = associations_p_l(1,i);
			state_landmark = id2state_landmark(id_landmark);
			temp_ass(:,end+1) = [state_pose;state_landmark];
			temp_meas(end+1) = Z_range(i);
		end
	end

	associations_p_l = temp_ass;
	Z_range = temp_meas;

	num_poses = length(state2id_pose);
	num_landmark = length(state2id_landmark);
	% poses in an array of 3x3 homogeneous transform matrices
	XR_guess=zeros(3,3,num_poses);

	for (pose_num=1:num_poses)
	    xr=initial_guess_poses(pose_num,:);
	    Xr=v2t(xr);
	    XR_guess(:,:,pose_num)=Xr;
	endfor;

	%now i need to find the initial guess of the landmarks using multilateration

	XL_guess = zeros(2,num_landmark);
	ranges = [];
	positions = [];


	for (lan_num=1:num_landmark)
		[rs,ps] = find_ranges(lan_num,associations_p_l,Z_range,XR_guess);
		land_pos = triangulate(ps,rs);
		XL_guess(:,lan_num) = land_pos;
	end

end


function [XR_ground_truth,XL_ground_truth,Zij_ground_truth,id2state_landmark_ground_truth] = load_ground_truth(file)
	XL_ground_truth = [];
	XR_ground_truth = [];
	id2state_landmark_ground_truth = [];
	Zij_ground_truth = zeros(3,3,1);
	flag = 0;
	f = fopen(file);
	l = fgetl(file);
	while ischar(l)
		d = strsplit(l);
		header = d{1};
		%get all initial guesses for the robot poses
		if strcmp(header,"VERTEX_XY") == 1
			id_landmark = str2double(d{2});
			x_l = str2double(d{3});
			y_l = str2double(d{4});
			XL_ground_truth(:,end+1) = [x_l;y_l];
			id2state_landmark_ground_truth(id_landmark) = size(XL_ground_truth,2);
		end
		if strcmp(header,"VERTEX_SE2") == 1
			x_r = str2double(d{3});
			y_r = str2double(d{4});
			theta_r = str2double(d{5});
			if flag == 0
				XR_ground_truth(:,:,1) = v2t([x_r;y_r;theta_r]);
				flag = 1;
			else	
				XR_ground_truth(:,:,end+1) = v2t([x_r;y_r;theta_r]);
			end
		end
		l = fgetl(file);
	end

	for(i=2:size(XR_ground_truth,3))
		XR_i = XR_ground_truth(:,:,i-1);
		XR_j = XR_ground_truth(:,:,i);
		dX = inv(XR_i)*XR_j;
		Zij_ground_truth(:,:,i-1) = dX;
	end
end

%function to flatten an homogenous transformation matrix matrix
% X: input homogenous transformation matrix matrix containing R = (r1 r2) and t
% flat: output flattened matrix (r1' r2' t')'

function flat = flatten(X)
	r1 = X(1:2,1);
	r2 = X(1:2,2);
	t = X(1:2,3);
	flat = [r1;r2;t];
end

% function to compute Multilateration solving a minimization problem
% input:
%	ps: (2xnum_poses) vector of poses from which the landmark has been observed
%	rs: (1xnum_poses) vector of ranges measurement made by poses i in ps with the correct order
% output
%	l: a guess for the landmark position which solves the minimization on the linear system Ax = b

function l = triangulate(ps,rs)
	num_points = length(rs);
	last_x = ps(1,num_points);
	last_y = ps(2,num_points);
	r_last = rs(num_points);
	A = zeros(num_points-1,2);
	b = zeros(num_points-1,1);
	for(i=1:num_points-1)
		xi = ps(1,i);
		yi = ps(2,i);
		ri = rs(i);
		A(i,1) = xi - last_x;
		A(i,2) = yi - last_y;
		b(i) = xi^2-last_x^2 + yi^2-last_y^2 + r_last^2-ri^2;
	end
	b = 0.5*b;
	damp_pinv = inv(A'*A+eye(2)*0.001)*A';

	l = damp_pinv*b;
end

% function that finds all the ranges measurement associated to a landmark
% input:
%	land_id: id of the target landmark
% 	ass: 2xnum_measurements. 
%      	 ass(:,k)=[p_idx,l_idx]' means the kth measurement
%        refers to an observation made from pose p_idx, that
%        observed landmark l_idx
% 	r_meas: vector of all the range measurements
% 	XR_guess: (3,3,num_poses) initial guesses of the robot poses
% output:
%	rs: vetcor of range measurements associated to land_id
%	px: vector of the poses from which the measurements in rs were taken

function [rs,px] = find_ranges(land_id,ass,r_meas,XR_guess)
	rs = [];
	px = [];
	for(i=1:length(r_meas))
		l_id = ass(2,i);
		if l_id == land_id
			z = r_meas(i);
			rs(end+1) = z;
			pose_id = ass(1,i);
			Xr = XR_guess(:,:,pose_id);
			euclidian_state = t2v(Xr);
			px(:,end+1) = euclidian_state(1:2);
		end
	end
end