function y = force1(c1,c1p,c2,c2p,c3,c3p,c4,c4p,c5,c5p,c6,c6p, ...
                    c7,c7p,c8,c8p,c9,c9p,c10,c10p,c11,c11p)

% force1.m
% This function receives receives the coordinates of the eight arm
% insertion points in the order of the ArmIP matrix of pointy.m, and 3
% forearm insertion points in the order of the FarmIP matrix. Along with 
% these, it receives the derivatives of each one of these points. Using
% this it produces velocity errors for each of the 10 muscles, and the 3D
% vectors along which each one of the muscles exerts force.
% The forces produced will be aimed at placing the hand at the specified
% target coordinates using the 2-level PCT algorithm.
% Notice there are 10 muscles, but 11 arm/forearm insertion points. From
% the 11 insertion points 1 corresponds to 2 muscles, and 2 muscles have
% both insertion points in the arm/forearm.
% This program is used by Arm10a.mdl.

% The coordinates are received in cm.
% The derivatives of the coordinates are in cm/s.
% The errors should be returned in cm/s.
% The ouput is a (10*3)+10 = 40 element column matrix. The first 30
% elements contain the unit action vectors for the muscles, pointing
% towards the shoulder (if shoulder to arm muscles) or the arm (if arm to
% forearm muscles). The remaining 10 elements are the velocity errors.

global ShouldIP Tlengths xgain vgain tension

% 1) Calculate current muscle lengths and the action vectors
% remember that ArmIP(5,:) = ArmIP(7,:)
RotArmIP = [c1; c2; c3; c4; c5; c6; c5; c7; c8];
RotFarmIP = [c9; c10; c11];
DtArmIP = [c1p; c2p; c3p; c4p; c5p; c6p; c5p; c7p; c8p];
DtFarmIP = [c9p; c10p; c11p];

v = zeros(10,3); % each row is a unit vector in the direction of action
v(1:7,:) = ShouldIP(1:7,:) - RotArmIP(1:7,:); % first 7 muscles of (41)
v(8,:) = ShouldIP(8,:) - RotFarmIP(1,:);  % biceps
v(9:10,:) = RotArmIP(8:9,:) - RotFarmIP(2:3,:); % brachialis and triceps

lengths = sqrt(sum(v.*v,2));

% 2) Obtain the derivative of muscle lengths (based on the back of (47))
%Dtlenghts = zeros(10,1);
Dtlengths(1:7,1) = -sum(v(1:7,:).*DtArmIP(1:7,:),2)./lengths(1:7); % first 7
Dtlengths(8,1) = -sum(v(8,:).*DtFarmIP(1,:),2)/lengths(8); % biceps
Dtlengths(9:10,1) = sum(v(9:10,:).*(DtArmIP(8:9,:)-DtFarmIP(2:3,:)),2)./lengths(9:10);

v = bsxfun(@rdivide,v,lengths+(1e-5)); % normalize the action vectors

% 3) Obtain the errors
%forces = zeros(3,11);
% a positive x_error means we need a negative velocity (muscle
% contraction), which will be effected if v_error > 0
x_error = xgain*(lengths - (Tlengths-tension));
v_error = vgain*max(x_error + Dtlengths,0);


% LEGACY CODE (from Arm9d):
% shoulder mounted muscles
% forces(:,1:8) = vgain*bsxfun(@times,v_error(1:8),v(1:8,:))';
% 
% % the following is because I made the crappy decision of considering
% % points 5 and 7 as the same point since they have the same location
% forces(:,5) = forces(:,5) + forces(:,7);
% forces(:,9) = forces(:,8); % forces(:,7) = biceps force on forearm
% 
% % next the "arm to forearm" muscles (brachialis and triceps)
% forces(:,7) = -v_error(9)*v(9,:)'; % brachialis arm
% forces(:,10) = v_error(9)*v(9,:)'; % brachialis farm
% forces(:,8) = -v_error(10)*v(10,:)'; % triceps arm
% forces(:,11) = v_error(10)*v(10,:)'; % triceps farm
% %forces = zeros(3,11);
% 
% % 3) Add the forces that avoid overstretch of muscles
% %stretch = 10000*exp(-40*abs(maxlen - lengths));
% % stretch = min(10000,(repmat(5,10,1)./abs(maxlen - lengths)).^4);
% % forces(:,1:6) = forces(:,1:6) + bsxfun(@times,stretch(1:6),v(1:6,:))';
% % forces(:,5) = forces(:,5) + stretch(7)*v(7,:)'; % the 7th muscle on the 5th IP
% % forces(:,7) = forces(:,7) + stretch(8)*v(8,:)'; % biceps on forearm IP 1
% % forces(:,8) = forces(:,8) - stretch(9)*v(9,:)'; % brachialis on arm
% % forces(:,9) = forces(:,9) + stretch(9)*v(9,:)'; % brachialis on farm
% % forces(:,10) = forces(:,10) - stretch(10)*v(10,:)'; % triceps on arm
% % forces(:,11) = forces(:,11) + stretch(10)*v(10,:)'; % triceps on farm
% 

y = [reshape(v,30,1); v_error];

