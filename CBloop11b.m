function y = CBloop11b(u)

% y = CBloop11b(u) implements the cerebellar main loop as described in
% scratch (38). This program is an adaptation of M1RBF2level.m, which
% refines the procedures of CBloop3.m in Arm7.
% The difference with CBloop11 is that there are several separate
% "cerebellar modules," each with its own particular timing for learning
% corrections (see (58) and (59)).

% For use with Arm11b.mdl

% The inputs to CBloop11 are:
% u(1) = simulation time  [s]
% u(2:12) = velocity errors (control signal) [cm/s]
% u(13:16) = shoulder quaternion (in the order of the shoulder sensor)
% u(17:20) = derivative of shoulder quaternion
% u(21) = elbow angle [rad]
% u(22) = elbow angular velocity [rad/s]
% u(23) = error (distance from hand to target) [cm]
% u(24) = positive part of error derivative [cm/s]


% From targen as global variables:
% shoulderQ = desired shoulder orientation [quat] (in 'ZXZ' order?)
% delta = desired elbow angle [rad]

% gammaC is still constant
% The indices for W are the TRANSPOSED of those in (29) and (30)


global MAXNF NF F W TLCS gammaC S epC ...
       states Times MAXNst Fis delta shoulderQ N3 N7 phi angF pc
       % avgSt T tauC CSsum M m % not currently used

% feature vector
v = [shoulderQ delta u(2:23)'];

% positive part of error derivative
CS = max(u(24),0);

step = u(1) - Times(MAXNst);
if step > 0  % because the function may be accessed many times
    % before the time advances
    
    % storing the new state
    states = circshift(states,[-1,0]);
    states(end,:) = v;
    Times = circshift(Times,-1);
    Times(MAXNst) = u(1);
    
    for i=1:N3+N7
        if S(i) == 0   % no unhandled spikes from previous steps
            % calculating probability to spike
            Pcs = pc*(cos(angF(i)*(u(1) - phi(i)))+1) / ...
                ( (1+exp(5-u(23)))*(1+exp(40-10*CS)) ); % see (59)
            
            if (Pcs > rand) && (u(1)-TLCS(i) > 0.2) % complex spike generated
                TLCS(i) = u(1);
                S(i) = 1;
                
                if NF(i) < MAXNF
                    NF(i) = NF(i) + 1;
                    %vDtauC = interp1(Times,states,max(Times(1),(u(1)-tauC)));
                    %F(NF,:) = vDtauC; % the new feature is the delayed input
                    disp(['New feature ' num2str(NF(i)) ' at time ' ...
                        num2str(u(1)) ' for module ' num2str(i)])
                end
            end
        else    % there are unhandled spikes from previous steps
            
            if CS <= 0 || (u(1) - TLCS(i)) > 0.35 % time to store the feature
                % obtaining the average control signal
                controls = states(Times>TLCS(i),6:16);
                meanC = mean(controls);
                
                base = (i-1)*MAXNF;
                if NF(i) < MAXNF
                    W(base+NF(i),:) = epC*meanC;
                    F(base+NF(i),:) = ...
                        interp1(Times,states,max(Times(1),1.5*TLCS(i)-0.5*u(1)));
                    % 1.5*TLCS-0.5*u(1) = TLCS - 0.5*(t - TLCS)
                    disp(['Feature ' num2str(base+NF(i)) ' received weights: ' ...
                        num2str(W(base+NF(i),:))]);
                else
                    FDel = interp1(Times,states,max(Times(1),1.5*TLCS(i)-0.5*u(1)));
                    dists = sum((bsxfun(@minus,F(base+1:base+MAXNF,:),FDel).^2),2);
                    [~,in] = min(dists);
                    F(base+in,:) = F(base+in,:) + 0.3*(FDel - F(base+in,:));
                    W(base+in,:) = W(base+in,:) + 0.3*(epC*meanC - W(base+in,:));
                    
                    disp(['Feature ' num2str(base+in) ' got blended'])
                    disp(num2str(u(1) - TLCS(i)))
                end
                
                S(i) = 0;
                
            end
        end
    end
end

dists = sum((bsxfun(@minus,F,v).^2),2);
dists = sqrt(MAXNF)*dists/norm(dists);

% Gaussian kernel
%errorF = 1/(1 + exp(-0.3*(CSsum-4)));
%errorF = 1/(1 + exp(-0.3*(CS-3)));
%Fis = (W')*exp(-gammaC.*dists);
Fis = min(max((W'),-3),3)*exp(-gammaC.*dists);
Fis = min(max(Fis,-3),3);
%Fis = errorF*(Fis + .9*(nFis - Fis));

% Linear kernel 
%Fis = 0.5*(W')*dists;
%nFis = (W')*max(1-10*dists,0);

% Fis = Fis + .9*(nFis - Fis);

%Fis = Fis/(norm(Fis) + .2);
% THE RIGHT WAY TO NORMALIZE Fis MAY BE WITH: sum(exp(-gammaC.*dists))
y = Fis;
%y = -2.5*ones(4,1);
%y = 1;

