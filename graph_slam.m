source "./help/geometry_helpers_2d.m"

%(minimal) size of pose and landmarks
global pose_dim=3;
global landmark_dim=2;


% retrieves the index in the perturbation vector, that corresponds to
% a certain pose
% input:
%   pose_index:     the index of the pose for which we want to compute the
%                   index
%   num_poses:      number of pose variables in the state
%   num_landmarks:  number of pose variables in the state
% output:
%   v_idx: the index of the sub-vector corrsponding to 
%          pose_index, in the array of perturbations  (-1 if error)
function v_idx=poseMatrixIndex(pose_index, num_poses, num_landmarks)
  global pose_dim;
  global landmark_dim;

  if (pose_index>num_poses)
    v_idx=-1;
    return;
  endif;
  v_idx=1+(pose_index-1)*pose_dim;
endfunction;


% retrieves the index in the perturbation vector, that corresponds to
% a certain landmark
% input:
%   landmark_index:     the index of the landmark for which we want to compute the
%                   index
%   num_poses:      number of pose variables in the state
%   num_landmarks:  number of pose variables in the state
% output:
%   v_idx: the index of the perturnation corrsponding to the
%           landmark_index, in the array of perturbations
function v_idx=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks)
  global pose_dim;
  global landmark_dim;
  if (landmark_index>num_landmarks)
    v_idx=-1;
    return;
  endif;
  v_idx=1 + (num_poses)*pose_dim + (landmark_index-1) * landmark_dim;
endfunction;

% error and jacobian of a measured range of the couple pose-landmark
% input:
%   Xr: the robot pose (3x3 homogeneous matrix) wrt world
%   Xl: the landmark pose (2x1 vector, 2d pose in world frame)
%   z:  range between the robot and landmark
% output:
%   e: 1x1 is the difference between prediction and measurement
%   Jr: 1x3 derivative w.r.t a the error and a perturbation on the
%       pose
%   Jl: 1x2 derivative w.r.t a the error and a perturbation on the
%       landmark
function [e,Jr,Jl]=errorAndJacobian(Xr,Xl,z)
   R=Xr(1:2,1:2);
   t=Xr(1:2,3);
   p_hat = R'*(Xl -t);  % Xr^(-1)*Xl
   z_hat=norm(p_hat); %prediction
   e=z_hat-z;
   Jr = zeros(1,3);
   J_icp = zeros(2,3);
   J_icp(1:2,1:2) = -R';
   J_icp(1:2,3) = R'*[0 1;-1 0]*Xl;
   Jr = (1/norm(p_hat))*p_hat'* J_icp;
   Jl= (1/norm(p_hat))*p_hat'*R';
endfunction;

% error and jacobian of a pose pose computed using the Chordal distance
% input:
%   Xri: the i-th robot pose (3x3 homogeneous matrix) wrt world
%   Xrj: the j-th pose (3x3 homogeneous matrix) wrt world
%   Zij:  Xij (3x3 homogenous matrix) displacement between Xj and Xi build through the odometry measurements
% output:
%   e: 3x1 is the difference between prediction and measurement
%   Jri: 6x3 derivative w.r.t a the error and a perturbation on the
%       pose Xi
%   Jl: 6x3 derivative w.r.t a the error and a perturbation on the
%       pose Xj
function [eij,Jri,Jrj]=errorAndJacobianPosePose(Xri,Xrj,Zij)
   Ri=Xri(1:2,1:2);
   ti=Xri(1:2,3);
   Rj=Xrj(1:2,1:2);
   tj=Xrj(1:2,3);
   % compute the relative displacemente between Xj and Xi
   Zij_hat = inv(Xri)*Xrj;
   % Since I am using the Chordal Distance I need to flatten the measurement and the prediction
   Zij_hat_flat = flatten(Zij_hat);
   Zij_flat = flatten(Zij);
   eij=Zij_hat_flat - Zij_flat;
   Jrj = zeros(6,3);
   % I compute the derivative wrt to each component in the perturbation vector associated to the pose Xj and build the final jacobian
   dh_dtx_j = [zeros(2,2),Ri'*[1;0]];
   dh_dty_j = [zeros(2,2),Ri'*[0;1]];
   dR_0 = [0 -1;1 0];
   dh_dtheta_j = [Ri'*dR_0*Rj Ri'*dR_0*tj];
   % Now I need to flatten each partial derivative in order to build the correct Jacobian for the chordal distance error
   Jrj(:,1) = flatten(dh_dtx_j);
   Jrj(:,2) = flatten(dh_dty_j);
   Jrj(:,3) = flatten(dh_dtheta_j);
   Jri = -Jrj;
endfunction;

% Function to compute the gamma term when using the the Geman-McLaure cost function
% input:
%   u: norm of the error
% output:
%   gamma_ : the gamma factor that needs to be absorbed in the information matrix 

function gamma_ = compute_gamma(u)
  drho_du = u/((1+u)^2);
  gamma_ = (1/u)*drho_du; 
end


% implementation of the boxplus
% applies a perturbation to a set of landmarks and robot poses
% input:
%   XR: the robot poses (3x3xnum_poses: array of homogeneous matrices)
%   XL: the landmark pose (2xnum_landmarks matrix of landmarks)
%   num_poses: number of poses in XR (added for consistency)
%   num_landmarks: number of landmarks in XL (added for consistency)
%   dx: the perturbation vector of appropriate dimensions
%       the poses come first, then the landmarks
% output:
%   XR: the robot poses obtained by applying the perturbation
%   XL: the landmarks obtained by applying the perturbation
function [XR, XL]=boxPlus(XR, XL, num_poses, num_landmarks, dx)
  global pose_dim;
  global landmark_dim;
  for(pose_index=1:num_poses)
    pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
    dxr=dx(pose_matrix_index:pose_matrix_index+pose_dim-1);
    XR(:,:,pose_index)=v2t(dxr)*XR(:,:,pose_index);
  endfor;
  for(landmark_index=1:num_landmarks)
    landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);
    dxl=dx(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,:);
    XL(:,landmark_index)+=dxl;
  endfor;
endfunction;

% implementation of the optimization loop with robust kernel
% applies a perturbation to a set of landmarks and robot poses
% input:
%   XR: the initial robot poses (3x3xnum_poses: array of homogeneous matrices) wrt world
%   XL: the initial landmark estimates (2xnum_landmarks matrix of landmarks)
%   Z:  the range measurements (1xnum_measurements)
%   Zij: the relative transformation between pose j and pose i (3x3xnum_poses)
%   associations_p_l: 2xnum_measurements. 
%                     associations_p_l(:,k)=[p_idx,l_idx]' means the kth measurement
%                     refers to an observation made from pose p_idx, that
%                     observed landmark l_idx
%   associations_p_p: 2xnum_measurements. 
%                     associations_p_p(:,k)=[pi_idx,pj_idx]' means the kth measurement
%                     refers to a relative trasformation between pose j and pose i
%   num_poses: number of poses in XR (added for consistency)
%   num_landmarks: number of landmarks in XL (added for consistency)
%   num_iterations: the number of iterations of least squares
%   damping:      damping factor (in case system not spd)
%   kernel_threshod: robust kernel threshold
%   cost_function : if 0 then simply normalizes square error when over a threshold, otherwise use the Geman-McLaure cost function in the minimizatioj

% output:
%   XR: the robot poses after optimization
%   XL: the landmarks after optimization
%   chi_stats: array 1:num_iterations, containing evolution of chi2
%   num_inliers: array 1:num_iterations, containing evolution of inliers
function [XR, XL, chi_stats,chi_stats_p_l,chi_stats_p_p, num_inliers]=doMultiICP(XR, XL, Z, Zij, 
							associations_p_l,
              associations_p_p, 
							num_poses, 
							num_landmarks, 
							num_iterations, 
							damping, 
							kernel_threshold,
              cost_function)
  global pose_dim;
  global landmark_dim;
  % Information Matrices for the Z and Zij measurement respectively
  sigma_range = 6;
  sigma_p_p = eye(6)*6;

  chi_stats=zeros(1,num_iterations);
  chi_stats_p_p = zeros(1,num_iterations);
  num_inliers=zeros(1,num_iterations);
  % size of the linear system
  system_size=pose_dim*num_poses+landmark_dim*num_landmarks; 
  for (iteration=1:num_iterations)
    H=zeros(system_size, system_size);
    b=zeros(system_size,1);
    chi_stats(iteration)=0;
    % First compute the error and accumulate in H and b all the terms due the range measurement error 
    for (measurement_num=1:size(Z,2))
      pose_index=associations_p_l(1,measurement_num);
      landmark_index=associations_p_l(2,measurement_num);
      z=Z(:,measurement_num);
      Xr=XR(:,:,pose_index);
      Xl=XL(:,landmark_index);
      [e,Jr,Jl] = errorAndJacobian(Xr, Xl, z);
      gamma_ = 1;
      % I need to handle differently the computation of some terms based on what cost function I am actually minimizing
      if cost_function == 0
        chi=e'*e;
        if (chi>kernel_threshold)
          e*=sqrt(kernel_threshold/chi);
          chi=kernel_threshold;
        else
          num_inliers(iteration)++;
        endif;
      else
        u = sqrt(e'*e);
        gamma_ = compute_gamma(u);
        chi = u^2/(2*(1+u^2));
      end
      chi_stats(iteration)+=chi;

      Hrr=Jr'*gamma_*sigma_range*Jr;
      Hrl=Jr'*gamma_*sigma_range*Jl;
      Hll=Jl'*gamma_*sigma_range*Jl;
      br=Jr'*gamma_*sigma_range*e;
      bl=Jl'*gamma_*sigma_range*e;

      % with the range measurement I am effecting the block in H wrt to: 
      %   pose_i x landmark_m
      %   pose_i x pose_i
      %   landmark_m x landmark_m
      %   landmark_m x pose_i
      % wrt to b I am effecting block:
      %   pose_i
      %   landmark_m

      pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
      landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);

      H(pose_matrix_index:pose_matrix_index+pose_dim-1,
	pose_matrix_index:pose_matrix_index+pose_dim-1)+=Hrr;

      H(pose_matrix_index:pose_matrix_index+pose_dim-1,
	landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Hrl;

      H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
	landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Hll;

      H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
	pose_matrix_index:pose_matrix_index+pose_dim-1)+=Hrl';

      b(pose_matrix_index:pose_matrix_index+pose_dim-1)+=br;
      b(landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=bl;

    endfor

    chi_stats_p_l(iteration) = chi_stats(iteration);

  %{  
    % Then compute the error and keep accumulating in H and b all the terms due the relative poses measurement error
    for (measurement_num=1:size(Zij,3))
      pose_index_i=associations_p_p(1,measurement_num);
      pose_index_j=associations_p_p(2,measurement_num);
      z=Zij(:,:,measurement_num);
      Xri=XR(:,:,pose_index_i);
      Xrj=XR(:,:,pose_index_j);
      [e,Jri,Jrj] = errorAndJacobianPosePose(Xri, Xrj, z);
      gamma_ = 1;
      if cost_function == 0
        chi=e'*e;
        if (chi>kernel_threshold)
          e*=sqrt(kernel_threshold/chi);
          chi=kernel_threshold;
        else
          num_inliers(iteration)++;
        endif;
      else
        u = sqrt(e'*e);
        gamma_ = compute_gamma(u);
        chi = u^2/(2*(1+u^2));
      end
      chi_stats(iteration)+=chi;
      chi_stats_p_p(iteration)+=chi;

      Hii=Jri'*gamma_*sigma_p_p*Jri;
      Hij=Jri'*gamma_*sigma_p_p*Jrj;
      Hjj=Jrj'*gamma_*sigma_p_p*Jrj;
      bri=Jri'*gamma_*sigma_p_p*e;
      brj=Jrj'*gamma_*sigma_p_p*e;

      % with the range measurement I am effecting the block in H wrt to: 
      %   pose_i x pose_j
      %   pose_i x pose_i
      %   pose_j x pose_j
      %   pose_j x pose_i
      % wrt to b I am effecting block:
      %   pose_i
      %   pose_j

      pose_matrix_index_i=poseMatrixIndex(pose_index_i, num_poses, num_landmarks);
      pose_matrix_index_j=poseMatrixIndex(pose_index_j, num_poses, num_landmarks);

      H(pose_matrix_index_i:pose_matrix_index_i+pose_dim-1,
  pose_matrix_index_i:pose_matrix_index_i+pose_dim-1)+=Hii;

      H(pose_matrix_index_i:pose_matrix_index_i+pose_dim-1,
  pose_matrix_index_j:pose_matrix_index_j+pose_dim-1)+=Hij;

      H(pose_matrix_index_j:pose_matrix_index_j+pose_dim-1,
  pose_matrix_index_j:pose_matrix_index_j+pose_dim-1)+=Hjj;

      H(pose_matrix_index_j:pose_matrix_index_j+pose_dim-1,
  pose_matrix_index_i:pose_matrix_index_i+pose_dim-1)+=Hij';

      b(pose_matrix_index_i:pose_matrix_index_i+pose_dim-1)+=bri;
      b(pose_matrix_index_j:pose_matrix_index_j+pose_dim-1)+=brj;

    endfor
    %}
    

    H+=eye(system_size)*damping;
    dx=zeros(system_size,1);

    % we solve the linear system, blocking the first pose
    % this corresponds to "remove" from H and b the locks
    % of the 1st pose, while solving the system

    dx(pose_dim+1:end)=-(H(pose_dim+1:end,pose_dim+1:end)\b(pose_dim+1:end,1));
    [XR, XL]=boxPlus(XR,XL,num_poses, num_landmarks, dx);
  endfor
endfunction


% plot landmarks and poses
%
%
%
function i = plotState(XL, XL_guess, XL_gt)
%plot landmarks
hold on;
plot(XL(1,:),XL(2,:),'b*',"linewidth",2);
hold on;
plot(XL_guess(1,:),XL_guess(2,:),'ro',"linewidth",2);
hold on;
plot(XL_gt(1,:),XL_gt(2,:),'g*',"linewidth",2);
hold on;
legend("estimate","initial guess","ground truth")
i = 1;
endfunction

function i = plotStateRobot(XR, XR_guess, XR_gt)
%plot landmarks
hold on;
plot(XR(1,3,:),XR(2,3,:),'b*',"linewidth",2);
hold on;
plot(XR_guess(1,3,:),XR_guess(2,3,:),'ro',"linewidth",2);
hold on;
plot(XR_gt(1,3,:),XR_gt(2,3,:),'g*',"linewidth",2);
hold on;
legend("estimate","initial guess","ground truth")
i = 1;
endfunction
