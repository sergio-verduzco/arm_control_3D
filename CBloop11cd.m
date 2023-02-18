function y = CBloop11cd(u)

% y = CBloop11cd(u) implements a variation of the cerebellar loop, using
% a combination of the error signals in CBloop11c and CBloop1d.
% In this version of CBloop errors are generated whenever the visual
% error (distance between the hand and the target) is increasng, but
% corrections are applied only to those muscles that were extending when
% they should have been contracting. Basically, the only difference with
% CBloop11d is that I uncommented  "(CS>0)" in line 81.
% The correction associated with that error is always an anticipative 
% strengthening of its contraction.

% As in CBloop11(c|d), we assume a long visual delay, so
% the error derivative is delayed by the value of the global variable vd.
% All the other inputs (except for the simulation time) are delayed by the
% value of d.

% For use with Arm11d.slx and params11d.m

% The inputs to CBloop11cd are:
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


global MAXNF NF F W TLCS gammaC S epC Tlengths maxlen N3 N7...
    states Times MAXNst Fis phi angF pc delta shoulderQ vd d 
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
       
    for i=1:N3+N7
        if S(i) == 0   % no unhandled spikes from previous steps
            % calculating probability to spike
            Pcs = pc*(cos(angF(i)*(u(1) - phi(i)))+1) / ...
                ( (1+exp(5-u(34)))*(1+exp(30-15*CS)) ); % see (59)
            
            if (Pcs > rand) && (u(1)-TLCS(i) > 0.2) % complex spike generated
                TLCS(i) = u(1);
                S(i) = 1;
                
                if NF(i) < MAXNF
                    NF(i) = NF(i) + 1;                 
                    disp(['New feature ' num2str(NF(i)) ' at time ' ...
                        num2str(u(1)) ' for module ' num2str(i)])
                end
            end
        else    % there are unhandled spikes from previous steps
            
            if CS <= 0  % time to store the feature
                %controls = zeros(1,11); controls(i) = 1; %correction involves ony one muscle
                
                % Finding NpLens past lengths used to generate the correction
                NpLens = 8;
                prevTimes = linspace(TLCS(i)-vd+d,TLCS(i)+0.7*(u(1)-TLCS(i))-vd+d,NpLens);
                prevLens = interp1(Times,states(:,28:38),prevTimes);
                % Approximating the mean derivative of the past lengths
                DpLs = bsxfun(@rdivide,prevLens(2:NpLens,:) - prevLens(1:NpLens-1,:),(prevTimes(2:NpLens)-prevTimes(1:NpLens-1)+1e-5)');
                DpL = mean(DpLs);
                % You generate a correction if length was increasing
                % when it should have been decreasing
                controls = max(mean(prevLens)-Tlengths',0).*max(DpL,0);

                base = (i-1)*MAXNF;
                if NF(i) < MAXNF
                    W(base+NF(i),:) = epC*controls;
                    timD = min(0.5*(u(1)-TLCS(i)),0.2);
                    timeC = max(Times(1),TLCS(i)-timD-vd+d);
                    F(base+NF(i),:) = ...
                        interp1(Times,states(:,1:27),timeC);
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
                
            else
                if ((u(1) - TLCS(i)) > 0.4 && CS > 0.1)|| u(34)<2 % this took too long; not storing
                    S(i) = 0;
                    NF(i) = NF(i)-1;
                    disp(['New feature ' num2str(NF(i)) ' for module' ...
                        num2str(i) ' got cancelled ']);
                end
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

