% params11cd.m
% Initialize the parameters of the model in Arm11d.slx when using 
% CBloop11cd. Run before each simulation (otherwise e_low won't update).

% Time is in seconds, length in centimeters, weight in kilograms

global first_sim  % could be set in reach_runner

% ----- these lines are not needed if calling from reach_runner
%first_sim = true;
% if first_sim
%     clear all; 
%     first_sim = true;
% end
% close all;
first_sim = false;
% -----

global ShouldIP maxlen xgain vgain p1 p2 tension last_t e_low ...
       LPfactor initL Larm Lfarm ArmIP FarmIP ...
       Tlengths MAXNF NF F W ...
       TLCS S epC gammaC states Times Nst MAXNst avgSt Fis ...
       shoulderQ vd d N3 N7 phi angF pc

%**************** GEOMETRIC PARAMETERS (same as piecePointy.m) **********
Larm = 30; % length of arm in cm
Rarm = 5; % radius of arm in cm
Lfarm = 26; % length of forearm in cm
Rfarm = 4; % radius of forearm in cm
Lhand = 7; % length of hand segment in cm
Rhand = 4; % radius of hand segment in cm

% Here are all the muscle insertion points in the trunk and shoulder, as
% labeled in (46) starting with point e
ShouldIP = [-5      0       0;     -11       0      -7; ...    % e,g
            -8     1        0;      -8      -1      0; ...      % i,k
            -14     6       1;      -10     -5      1; ...      % m,o
            -14     6       -12;    -10     -5      -12; ...    % q,s
            -2      2      1];     % t


% Here are all the muscle insertion points in the arm, as labeled in (46),
% starting with point d
ArmIP = [2      0       -5;     -2      0       -5; ...     % d,f
         0      2       0;      0       -2      0; ...      % h,j
         -1     1       -6;     -1      -1      -6; ...     % l,n
         -1     1       -5;     -1      -1      -5; ...     % p,r
         0      -1.5    -6;      0      1.5     -15];       % v,x
     
% Here are the muscle insertion points in the forearm, corresponding to the
% biceps, the triceps, and the brachialis.
FarmIP = [0     1.5     -Larm-5;      0       -0.5    -Larm+3; ...
          -1    1.5     -Larm-4];
     
% Here are the initial and final points of the restricted bending lines for
% each one of the muscles 1-8. The values for muscles 3,4 are not used
p1 = [3     -2      2;      1       3       -2; ...
      0     2       1;      0       -2      1; ...
      -2    3       2;     -1       -5      -7; ...
      -2    3       2;     -1       -5      -7];
  
p2 = [3     2       2;      1       -3      -2; ...
      0     2       -1;     0       -2      -1; ...
      -3    3       -6;     0       -4      2; ...
      -3    3       -6;     0       -4      2];
  
maxlen = 0.2+[16.8;19.2;10;10;21.4;20.2;24.2;21.9;36.2;26.5;11]; % maximum muscle lengths
 
 %********************************************

%############# MASS PROPERTIES ##############
dens = 0.001; % density of the arm in kg/cm^3 (density of water)
Marm = dens*pi*(Rarm^2)*Larm; % mass of the arm in kg
Mfarm = dens*pi*(Rfarm^2)*Lfarm; % mass of the forearm in kg
Mhand = dens*pi*(Rhand^2)*Lhand; % mass of hand in kg

% inertia tensor for arm:
Iarm = diag([(Marm/4)*(Rarm^2 + (Larm^2)/3), ...
       (Marm/4)*(Rarm^2 + (Larm^2)/3), (Marm/2)*(Rarm^2)]);
Iarm = Iarm*10^(-4); % transforming from kg*cm^2 to kg*m^2
% inertia tensor for forearm:
Ifarm = diag([(Mfarm/4)*(Rfarm^2 + (Lfarm^2)/3), ...
        (Mfarm/4)*(Rfarm^2 + (Lfarm^2)/3), (Mfarm/2)*(Rfarm^2)]);
Ifarm = Ifarm*10^(-4); % transforming from kg*cm^2 to kg*m^2
% inertia tensor for hand:
Ihand = diag([(Mhand/4)*(Rhand^2 + (Lhand^2)/3), ...
        (Mhand/4)*(Rhand^2 + (Lhand^2)/3), (Mhand/2)*(Rhand^2)]);
Ihand = Ihand*10^(-4); % transforming from kg*cm^2 to kg*m^2

%#############################################

%()()()()()()()()()() OTHER PARAMETERS ()()()()()()()()()()
xgain = 2; %.9; % xgain and vgain are for forces6.m
vgain = 1; %.2;
tension = .05*maxlen/max(maxlen);    % decreased from the target lengths
vd = 0.15; % visual delay
%()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
 
%*********** MUSCLE PARAMETERS *************
% Eq. 3
mu=0.06;  % dependence of muscle's threshold length on velocity (adimensional)
%mu=0;  % debugging value
d=0.025;  % reflex delay (in s)
%d=0;  % debugging reflex delay 

% Eq. 4
rhoP=14.9; % N (Newtons)
% rhoD=14.9; % N (Newtons)
% rhoB=13.1; % N (Newtons)
% rhoT=18.8; % N (Newtons)
         % rho the muscle's physiological cross-sectional area scaled by
         % the factor 1 (N/cm^2). The cross sectional areas for the arm
         % muscles are:
         % biceps short head = 2.1 cm^2; biceps long head = 11 cm^2;
         % deltoid = 14.9 cm^2; pectoralis = 14.9 cm^2;
         % triceps lateral head = 12.1 cm^2;
         % triceps long head = 6.7 cm^2 .
c=1.12; % cm^-1

% Eq. 5
tau = 0.01; % s
%tau = .1; % debugging value
tausq = tau^2;  % s^2

% Eq. 6
f1 = 0.82; % (adimensional)
f2 = 0.5;  % (adimensional)
f3 = 0.43;  % (adimensional)
f4 = 0.58;  % s/cm
           
kB = 1.0;  % N/cm
kD = 2.58;
kP = 2.58;
kT = 1.5;
           % The passive stiffness of the muscle is assumed to vary
           % linearly with physiological cross-sectional area. The scaled
           % values of k are:
           % biceps short head = .365 N/cm; biceps long head = 1.909 N/cm;
           % deltoid = 2.58 N/cm; pectoralis = 2.58 N/cm
           % triceps lateral head = 2.09 N/cm;
           % triceps long head = 1.163 N/cm

% these resting lengths are what pointy calculates given the angles
% alpha=-pi/4, beta=pi/4, gamma=0, delta=pi/4.
rlen = [13.8992; 12.5757; 9.1104; 7.1414; 17.4868; 17.7487; 18.8370; ...
        17.8958; 33.5286; 22.2536; 16.9666];
 %********************************************
 
%0000000000000000000000 L0W-PASS FILTER PARAMETERS 000000000000000000000
last_t = 0;  % time when e_low was last updated
e_low = zeros(11,1);  % low-pass filtered error signal for all muscles
LPtau = .2;   % time constant of the filter [s]
LPfactor = (1-exp(-1))/LPtau;  % factor for numerical integration of e_low
%00000000000000000000000000000000000000000000000000000000000000000000000

%###################% GENERATING VELOCITY TRAJECTORY ####################
% these are the initial lengths, copied from a run of piecePointy
% with alpha=pi/300,beta=pi/300,gamma=0,delta=0.
initL = [15.3170;  8.5450;  8.0830;  8.0415;  15.8521;  14.8531; ...
         15.8538;  14.8197; 36.0393; 21.0238; 19.0263];
%########################################################################     
     
 Tlengths = zeros(11,1);
 shoulderQ = zeros(1,5);
 
 %_________________________________________________________
 %;;;;;;;;;;;;;;;;;;;;;;; CEREBELLUM ;;;;;;;;;;;;;;;;;;;;;; 
 
 N3 = 6;    % number of 3 Hz modules
 N7 = 3;    % number of 7 Hz modules
 MAXNF = 50; % maximum number of features per module
 InpD = 38; % dimensionality of the input vector (without time and IO)
 OutD = 11; % dimensionality of the output vector (number of DCN cells)
 
 if first_sim
     % THE NEXT 3 LINES ERASE LEARNED FEATURES
     NF = zeros(N3+N7,1); % number of features stored in F matrix for module i
     F = zeros(N3*N7*MAXNF,InpD-11); % feature matrix
     W = zeros(N3*N7*MAXNF,OutD); % feature weights matrix
 end
 
 TLCS = -0.1*ones(N3+N7,1); % onset time of the last complex spike for module i
 %TtauC = T+tauC;
 S = zeros(N3+N7,1); % "unhandled spike" banner
 epC = .2   ; % controls the speed of change in the W weights 
 gammaC = 20*ones(N3*N7*MAXNF,1); % controls the width of the feature kernels 
 MAXNst = 3000; % maximum number of previous states stored
 Nst = 1; % number of previous states stored 
 
 %init = [pi/4 pi/4 pi/4 pi/4 pi/2 pi pi pi pi 0 0 0 0 pi 0 0];
 init = [.9 .2 .2 .2  pi/2 .5*ones(1,11) 0 0 0 1 0 0 0 0 pi 0 10 initL'];
 if length(init) == InpD
    states = repmat(init,MAXNst,1); % to save previous input states
 else
     disp('Error: Wrong initial statesize at params11d.m');
 end
 avgSt = init'; % time average of the state
 Times = (-(MAXNst-1):1:0)'; % time at which states(i) was saved 
 Fis = zeros(11,1);
 phi3 = (0:N3-1)/(3*N3);    % temporal phase of 3 Hz modules
 phi7 = (0:N7-1)/(7*N7);    % temporal phase of 7 Hz modules
 phi = [phi3 phi7];
 angF = [repmat(6*pi,N3,1);repmat(14*pi,N7,1)];
 pc = 0.01; % adjusts the probability of complex spikes
 %;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 %_________________________________________________________ 
         
disp('params11cd.m loaded')

     