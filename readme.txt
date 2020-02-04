the multipoint files are adaptation of multipoints 3d registration in 2d and using as measurement the range.
therefore the jacobians are different
--- measurement is z = norm(Xr^(-1)*Xl) because Xr is robot wrw and measurement are taken from the robot so i first convert the point from world to robot and then take the norm of this vector which is the distance of the landmar

--- the h(x[+]deltaX) = norm(v2t(-deltaXr)*Xr^(-1)*(Xl+deltaXl))

p_hat = Xr^(-1)*Xl
dR_minus is the derivative of R(-deltaTheta) = (cos(deltaTheta) sin(deltaTheta); -sin(deltaTheta) cos(deltaTheta))

(d(h(.))/d(deltaTr))(0) = (1/norm(p_hat))*(p_hat)'*(-R')
(d(h(.))/d(deltaTheta))(0) = (1/norm(p_hat))*(p_hat)'*(R'*dR_minus(0)*Xl)
(d(h(.))/d(deltaXl))(0) = (1/norm(p_hat))*(p_hat)'*R'


TO LAUNCH THE PROGRAM SIMPLY RUN THE main.m script
