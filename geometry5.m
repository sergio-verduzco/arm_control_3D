function y = geometry5(c1,c1p,c2,c2p,c3,c3p,c4,c4p,c5,c5p,c6,c6p,c7,c7p, ...
             c8,c8p,c9,c9p,c10,c10p,c11,c11p,c12,c12p,c13,c13p)

% geometry5.m
% This function receives receives the coordinates of the ten arm
% insertion points in the order of the ArmIP matrix of piecePointy.m, and 3
% forearm insertion points in the order of the FarmIP matrix. Along with 
% these, it receives the derivatives of each one of these points. Using
% this it produces the muscle lengths, the muscle velocities, and control 
% signals for each of the 11 muscles, and the 3D vectors along which the 
% muscles exert force.

% The controls produced will be aimed at placing the hand at the specified
% target coordinates using the velocity trajectory method.
% Notice there are 11 muscles, but 13 arm/forearm insertion points. This is
% because 2 muscles have both insertion points in the arm/forearm.

% The coordinates are received in cm.
% The derivatives of the coordinates are in cm/s.
% The ouput is a 11+11+11+(11*3) = 66 element column matrix. The first 33
% elements contain the lengths, contraction speeds, and velocity errors of
% the muscles. The remaining 33 are the unit action vectors for the 
% muscles, pointing towards the shoulder (if shoulder to arm muscles) or 
% the arm (if arm to forearm muscles). 

% This program is used by Arm11.mdl.

global ShouldIP Tlengths xgain vgain vgain2 maxlen p1 p2 tension initL ...
       x0 SM

% 1) Calculate current muscle lengths and the action vectors
% remember that ArmIP(5,:) = ArmIP(7,:)
RotArmIP = [c1; c2; c3; c4; c5; c6; c7; c8; c9; c10];
RotFarmIP = [c11; c12; c13];
DtArmIP = [c1p; c2p; c3p; c4p; c5p; c6p; c7p; c8p; c9p; c10p];
DtFarmIP = [c11p; c12p; c13p];

bend = zeros(8,4);
lengths = zeros(11,1);
v = zeros(11,3); % each row is a unit vector in the direction of action
for i = 1:8
    % muscles 1,2,5,6,7,8 may be wrapped around a bending line
    bend(i,:) = piece5(ShouldIP(i,:),RotArmIP(i,:),p1(i,:),p2(i,:));
    
    if bend(i,1) > 0  && i~=3 && i~=4 % if muscle wraps around bending line
        lengths(i) = sqrt(sum((ShouldIP(i,:) - bend(i,2:4)).^2)) + ...
                    sqrt(sum((bend(i,2:4) - RotArmIP(i,:)).^2));
        v(i,:) = bend(i,2:4) - RotArmIP(i,:);
    else        
        lengths(i) = sqrt(sum((ShouldIP(i,:) - RotArmIP(i,:)).^2));
        v(i,:) = ShouldIP(i,:) - RotArmIP(i,:); 
    end
end

lengths(9) = sqrt(sum((ShouldIP(9,:) - RotFarmIP(1,:)).^2,2));
lengths(10:11) = sqrt(sum((RotArmIP(9:10,:) - RotFarmIP(2:3,:)).^2,2));

v(9,:) = ShouldIP(9,:) - RotFarmIP(1,:);  % biceps
v(10:11,:) = RotArmIP(9:10,:) - RotFarmIP(2:3,:); % brachialis and triceps

% 2) Obtain the derivative of muscle lengths (based on the back of (47))
%Dtlenghts = zeros(10,1);
Dtlengths(1:8,1) = -sum(v(1:8,:).*DtArmIP(1:8,:),2)./lengths(1:8); % first 8
Dtlengths(9,1) = -sum(v(9,:).*DtFarmIP(1,:),2)/lengths(9); % biceps
Dtlengths(10:11,1) = sum(v(10:11,:).*(DtArmIP(9:10,:)-DtFarmIP(2:3,:)),2)./lengths(10:11);

v = bsxfun(@rdivide,v,lengths+(1e-5)); % normalize the action vectors

% 3) Obtain the errors
E = lengths - Tlengths;
% dir = sum(Dtlengths.*E)/(norm(Dtlengths)*norm(E)+1e-5);
%x_error = xgain*(1.5-dir)*(E + tension);
x_error = xgain*(E+tension);
v_error = vgain*(Dtlengths + x_error);
%x_error = xgain*(lengths - (Tlengths-tension));
%v_error = vgain*(Dtlengths - vgain2.*(sin(pi*(lengths-initL)./x0).^3));
%v_error = vgain*(Dtlengths - (sin(pi*(lengths-initL)./x0).^3));
%v_error = vgain*Dtlengths;
% control = x_error + v_error;
%control = v_error;

%control = control + SM*control; % antagonist inhibition
%v_error = v_error + SM*v_error; % antagonist inhibition

y = [lengths; Dtlengths; v_error; reshape(v,33,1)];


% legacy code from force6.m
% forces = zeros(3,13);
% a positive x_error means we need a negative velocity (muscle
% contraction), which will be effected if v_error > 0
% x_error = xgain*(lengths - (Tlengths - tension));
% v_error = vgain*max(x_error + Dtlengths,0);

%v_error = 1*ones(11,1);

%forces(:,1) = muscle_force(v_error,length(1),speed(1))*v(:,1)'

% first we must specify forces to arm IPs
% forces(:,1:8) = bsxfun(@times,v_error(1:8),v(1:8,:))'; % first 8 of (46)
% forces(:,9) = -v_error(10)*v(10,:)';   % triceps arm IP
% forces(:,10) = -v_error(11)*v(11,:)';   % brachialis arm IP
% forces(:,11) = v_error(9)*v(9,:)';     % biceps farm IP
% forces(:,12) = v_error(10)*v(10,:)';   % triceps farm IP
% forces(:,13) = v_error(11)*v(11,:)';   % brachialis farm IP
%forces(:,:) = zeros(3,13);
%forces(:,[1:11 13]) = zeros(3,12)
%forces(:,11) = v(9,:)';

% %Add the forces that avoid overstretch of muscles
% stretch = 100*exp(-40*max(lengths-maxlen,zeros(11,1)));
% %stretch = min(1000,(repmat(2,11,1)./abs(maxlen - lengths)).^2);
% forces(:,1:8) = forces(:,1:8) + bsxfun(@times,stretch(1:8),v(1:8,:))';
% forces(:,9) = forces(:,9) - stretch(10)*v(10,:)'; % triceps on arm
% forces(:,10) = forces(:,10) - stretch(11)*v(11,:)'; % brachialis on arm
% forces(:,11) = forces(:,11) + stretch(9)*v(9,:)'; % biceps on farm
% forces(:,12) = forces(:,12) + stretch(10)*v(10,:)'; % triceps on farm
% forces(:,13) = forces(:,13) + stretch(11)*v(11,:)'; % triceps on farm



