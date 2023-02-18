function y = CBloop11(u)

% y = CBloop11(u) implements the cerebellar main loop as described in
% scratch (38). This program is an adaptation of M1RBF2level.m, which
% refines the procedures of CBloop3.m in Arm7.
% The difference with CBloop is that the weights are calculated using an
% average of the M1 lambdas instead of using the contraction forces.

% For use with Arm11a.mdl

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
       states Times MAXNst muscles Nresp Fis muscAct mini in ...
       delta shoulderQ
       % avgSt T tauC CSsum M m % not currently used

% rescaling and normalizing the inputs
v = [shoulderQ delta u(2:23)'];

CS = max(u(24),0);

step = u(1) - Times(MAXNst);
if step > 0  % because the function may be accessed many times
    % before the time advances
    
    % storing the new state
    states = circshift(states,[-1,0]);
    states(end,:) = v;
    Times = circshift(Times,-1);
    Times(MAXNst) = u(1);
    
    if S == 0   % no unhandled spikes from previous steps
        if CS > 4 && (u(1)-TLCS > 0.2) % if we just received a few spikes
            TLCS = u(1);
            %CSsum = CS;
            S = 1;
            muscAct = u(2:12)'; % u(2:12) are the velocity errors
            Nresp = 1;
            %[mini in] = min(dists);
            if NF < MAXNF %&& (mini > .2)
                NF = NF + 1;
                %vDtauC = interp1(Times,states,max(Times(1),(u(1)-tauC)));
                %F(NF,:) = vDtauC; % the new feature is the delayed input
                disp(['New feature ' num2str(NF) ' at time ' num2str(u(1))])
            end
        end
    else    % there are unhandled spikes from previous steps
        %XT = min(.2,T*(1+(CSsum/20)));
        %if ((u(1) - TLCS >= T) && (CS < 6)) ||  (CS==0)  % time for plasticity!
        if CS > 0
            muscAct = muscAct + u(2:12)';
            Nresp = Nresp + 1;
        else
            %muscles([1 2]) = (muscAct([1 2])/(muscAct(1)+muscAct(2)+1e-5)-.5);
            %muscles([3 4]) = (muscAct([3 4])/(muscAct(3)+muscAct(4)+1e-5)-.5);
            muscles = muscAct/Nresp;
            if NF < MAXNF %&& (mini > .2)
                %W(NF,:) = W(NF,:) + epC*muscles;
                %errorF = 1/(1 + exp(-.2*(CSsum - 10)));
                %W(NF,:) = W(NF,:) + epC*errorF*muscles;
                %gammaC(NF) = 0.05 + 4*errorF;
                W(NF,:) = epC*muscles;
                % W = Weights (Corrections)
                F(NF,:) = interp1(Times,states,max(Times(1),1.5*TLCS-0.5*u(1)));
                % F = Features (Triggers weights, or corrections)
                % 1.5*TLCS-0.5*u(1) = TLCS - 0.5*(t - TLCS)
                disp(['Feature ' num2str(NF) ' received weights: ' ...
                    num2str(W(NF,:))]);
                %gammaC(NF) = 4*(1.02 - errorF);
            else
                FDel = interp1(Times,states,max(Times(1),1.5*TLCS-0.5*u(1)));                
                dists = sum((bsxfun(@minus,F,FDel).^2),2);                
                [mini,in] = min(dists);
                F(in,:) = F(in,:) + 0.3*(FDel - F(in,:));
                W(in,:) = W(in,:) + 0.3*(epC*muscles - W(in,:));                
                %F(in,:) = F(in,:)/(norm(F(in,:)+1e-5));
                %errorF = 1/(1 + exp(-.2*(CSsum - 10)));
                %W(in,:) = W(in,:) + 0.3*(epC*errorF*muscles-W(in,:));
                %W(in,:) = W(in,:) + epC*errorF*(M-W(in,:)).*(W(in,:)-m).*muscles;
                %gammaC(in) = gammaC(in) + 0.3*(0.05 + 4*errorF - gammaC(in));
                disp(['Feature ' num2str(in) ' got blended'])
                disp(num2str(u(1) - TLCS))
            end
                    
            S = 0;
            %CSsum = 0;
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

