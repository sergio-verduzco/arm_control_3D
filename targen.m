function Y = targen(t)

% Y = targen(t) generates the (x,y,z) coordinates of the target for the
% reaching task. These coordinates may change over time. targen also
% generates the quaternion corresponding to the shoulder orientation when
% the hand is at the target, and passes it as a global variable.
% The coordinates produced are in cm.

% Input: Current simulation time (t)
% Output: Target coordinates

global Larm Lfarm ArmIP FarmIP ShouldIP delta gamma beta alpha ...
       Tlengths p1 p2 shoulderQ x0 initL vgain2

%>>>>>>>>>>>>>>>>>>>>>> Motor Cortex <<<<<<<<<<<<<<<<<<<<<<
% target coordinates (careful to avoid Gimble locks)
xt = 30 + 0*sin(1e-5*t);
yt = 30 + 0*cos(1e-5*t);
zt = -10 - 0*cos(1e-5*t);

Y = [xt,yt,zt];
D = norm(Y);

% From coordinates to angles as in (44)
delta = pi - acos((Larm^2 + Lfarm^2 - D^2)/(2*Larm*Lfarm));
gamma = 0;
beta = acos(-zt/D) - acos((D^2 + Larm^2 - Lfarm^2)/(2*D*Larm));
alpha = asin(-xt/sqrt(xt^2 + yt^2));

% Now we'll rotate the insertion points by the Euler angles in reverse
% order 
% This magic R matrix (from (42)) does it all in one step
sa = sin(alpha); ca = cos(alpha);
sb = sin(beta);  cb = cos(beta);
sc = sin(gamma); cc = cos(gamma);

R = [ca*cc-sa*cb*sc, -ca*sc-sa*cb*cc, sa*sb;
     sa*cc+ca*cb*sc, -sa*sc+ca*cb*cc, -ca*sb;
     sb*sc,           sb*cc,          cb];
 
RotArmIP = (R*ArmIP')';

% For the forearm we also need the delta rotation
RotFarmIP = bsxfun(@plus,[0 0 -Larm], ...
         ([1 0 0; 0 cos(delta) -sin(delta); 0 sin(delta) cos(delta)]* ...
          bsxfun(@minus,FarmIP,[0 0 -Larm])')');
RotFarmIP = (R*RotFarmIP')';

% calculate the muscle lengths
bend = zeros(8,4);
for i = 1:8
    % muscles 1,2,5,6,7,8 may be wrapped around a bending line
    if sum(i == [1 2 5 6 7 8]) > 0
        bend(i,:) = piece5(ShouldIP(i,:),RotArmIP(i,:),p1(i,:),p2(i,:));
    end
    
    if bend(i,1) > 0   % if muscle wraps around bending line
        Tlengths(i) = sqrt(sum((ShouldIP(i,:) - bend(i,2:4)).^2)) + ...
                    sqrt(sum((bend(i,2:4) - RotArmIP(i,:)).^2));
    else        
        Tlengths(i) = sqrt(sum((ShouldIP(i,:) - RotArmIP(i,:)).^2));
    end
end

Tlengths(9) = sqrt(sum((ShouldIP(9,:) - RotFarmIP(1,:)).^2,2));
Tlengths(10:11) = sqrt(sum((RotArmIP(9:10,:) - RotFarmIP(2:3,:)).^2,2));

%()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

%###################% GENERATING VELOCITY TRAJECTORY ####################
% these are the initial lengths, copied from a run of piecePointy
% with alpha=pi/300,beta=pi/300,gamma=0,delta=0.
x0 = Tlengths - initL;
% Ldiffs = (1/100)*x0*(0:100);
% Ldiffs = [bsxfun(@minus,Ldiffs,x0), Ldiffs, bsxfun(@plus,x0,Ldiffs)];
% fvals = sin(pi*bsxfun(@rdivide,Ldiffs,x0)); % all rows are the same ^_^
vgain2 = abs(x0)/2;  % peak of velocity trajectory
%########################################################################

%================  CALCULATING SHOULDER QUATERNION  =================
% ca = cos(alpha/2); sa = sin(alpha/2); 
% cb = cos(beta/2);  sb = sin(beta/2);
% cg = cos(gamma/2); sg = sin(gamma/2);
% % Wikipedia formula...
% v1 = ca*cb*cg + sa*sb*sg;
% v2 = sa*cb*cg - ca*sb*sg;
% v3 = ca*sb*cg + sa*cb*sg;
% w = ca*cb*sg - sa*sb*cg;
% %shoulderQ = [v1;v2;v3;w];
% shoulderQ = [t v2 v3 w v1];
shoulderQ = angle2quat(alpha,beta,gamma,'ZXZ');
%=====================================================================
