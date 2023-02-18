function y = CBloop11d(u)

% y = CBloop11d(u) implements a variation of the cerebellar loop, using the
% ideas of 7/8/13 in log.txt.
% In this version of CBloop errors are detected and corrected at the level
% of individual muscles. Whenever a muscle that is longer than its
% equilibrium value gets longer, an error is generated. The correction
% associated with that error is always an anticipative strengthening of its
% contraction.

% As in CBloop11c, we assume a long visual delay, so
% the error derivative is delayed by the value of the global variable vd.
% All the other inputs (except for the simulation time) are delayed by the
% value of d.

% For use with Arm11d.mdl

% The inputs to CBloop11d are:
% u(1) = simulation time  [s]
% u(2:12) = muscle lengths [cm]
% u(13:23) = velocity errors (control signal) [cm/s]
% u(24:27) = shoulder quaternion (in the order of the shoulder sensor)
% u(28:31) = derivative of shoulder quaternion
% u(32) = elbow angle [rad]
% u(33) = elbow angular velocity [rad/s]
% u(34) = delayed error (distance from hand to target) [cm]
% u(35) = delayed positive part of error derivative [cm/s]


% From targen as global variables:
% shoulderQ = desired shoulder orientation [quat] (in 'ZXZ' order?)
% delta = desired elbow angle [rad]

% gammaC is still constant
% The indices for W are the TRANSPOSED of those in (29) and (30)


global MAXNF NF F W TLCS gammaC S epC Tlengths len_rat...
    states Times MAXNst Fis delta shoulderQ vd d nLenAvg avgLen maxlen
% avgSt T tauC CSsum M m % not currently used

% feature vector
v = [shoulderQ delta u(13:34)'];

% positive part of error derivative
CS = max(u(35),0);

% muscle lengths
lengths = u(2:12)';

step = u(1) - Times(MAXNst);
if step > 0  % because the function may be accessed many times
    % before the time advances
    
    % storing the new state
    states = circshift(states,[-1,0]);
    states(end,:) = [v, lengths];
    Times = circshift(Times,-1);
    Times(MAXNst) = u(1);
    
    avgLen = mean(states(MAXNst-nLenAvg:MAXNst,28:38)); % past average length
    
    for i=1:11
        if S(i) == 0   % no unhandled spikes from previous steps
            % if the muscle is getting longer and it's supposed to contract
            % and we haven't recently generated an error for this muscle and
            % the error is increasing
            if (lengths(i) > avgLen(i)) && (Tlengths(i) < lengths(i)-.0) ...
                    && (u(1)-TLCS(i) > 0.2) && (CS > 0)
                len_rat(i) = (lengths(i) - Tlengths(i)); %/maxlen(i);
                TLCS(i) = u(1);
                S(i) = 1;
                if NF(i) < MAXNF
                    NF(i) = NF(i) + 1;
                    disp(['New feature ' num2str(NF(i)) ' at time ' ...
                        num2str(u(1)) ' for muscle ' num2str(i)])
                end
            end
        else    % there are unhandled spikes from previous steps
            if Tlengths(i) >= lengths(i) || (u(1) - TLCS(i)) > 0.3  % time to store the feature
                %if CS <= 0  % time to store the feature
                    controls = zeros(1,11); controls(i) = 1;
                    base = (i-1)*MAXNF;
                    if NF(i) < MAXNF
                        W(base+NF(i),:) = (epC*len_rat(i))*controls;
                        F(base+NF(i),:) = ...
                            interp1(Times,states(:,1:27),max(Times(1),TLCS(i)-0.5*(u(1)-TLCS(i))-vd+d));
                        disp(['Feature ' num2str(base+NF(i)) ' received weights: ' ...
                            num2str(W(base+NF(i),:))]);
                    else
                        FDel = interp1(Times,states(:,1:27),max(Times(1),1.5*TLCS(i)-0.5*u(1)-vd+d));
                        dists = sum((bsxfun(@minus,F(base+1:base+MAXNF,:),FDel).^2),2);
                        [~,in] = min(dists);
                        F(base+in,:) = F(base+in,:) + 0.3*(FDel - F(base+in,:));
                        W(base+in,:) = W(base+in,:) + 0.3*(epC*controls - W(base+in,:));
                        
                        disp(['Feature ' num2str(base+in) ' got blended'])
                        disp(num2str(u(1) - TLCS(i)))
                    end
                    
                    S(i) = 0;
                    
%                 else
%                     if (u(1) - TLCS(i)) > 0.4 % this took too long; not storing
%                         S(i) = 0;
%                         NF(i) = NF(i)-1;
%                         disp(['New feature ' num2str(NF(i)) ' for module' ...
%                             num2str(i) ' got cancelled ']);
%                     end
%                 end
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

