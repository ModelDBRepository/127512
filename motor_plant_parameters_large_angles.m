% (c) Written By Erez Simony 2010, code for the model described in:  
% Simony, E., Bagdasarian K, Herfst L., Brecht M., Ahissar E, Golomb D. 
% Temporal and spatial characteristics of vibrissa responses to motor commands (2010). 
% Journal of Neuroscience, In press.

function motor_plant_parameters_large_angles


global N vib_num resting_angles intrinsic_muscle_set force_factor MN_spikes_times
global Lf s w Mh Mf C I0
global D zeta_up zeta_dn zeta_gr Kup Kdn Kgr
global r0 tauc taur 
global flow fhigh fl_stat 
global t_start t_step t_end time_shift

%%%%%%%%%%%%%    [ Common] %%%%%%%%%%%%%%%%%%%%%%%%
fl_stat=1;                                % Force-Length (Curve): 1:with FL curve; 0:w/o FL curve
N=5;                                      % Number of whiskers in a row
vib_num=4;                                % vibrissa number - most posterior=1.
resting_angles=[75 75 75 75 75 ];         % Whisker Resting angles - posterior resting angle: resting_angles(1)
intrinsic_muscle_set=[1;1;1;1;1;1];           % muscle inervations: 1-active muscle, 0-not active
amp=1.0877;                                    % muscle force amplitude      
force_factor=amp*(ones(1,length(intrinsic_muscle_set)));   %equal force amplitudes [amp amp amp amp amp amp];
% force_factor=[1 1.2 1.3 1.1 ];                           %Not equal force amplitudes [amp1 amp2 amp3 amp4 amp5 amp6];
MN_spikes_times=ones(N+1,1)*[0 15 30];
% MN_spikes_times=ones(N+1,1)*[0:12:100];

% %%% Bio-Mechanical Parametsrs (Pad)   [ Common]
Lf=4;                                     % follicle length (mm)
s=2;                                      % distance between 2 follicles (mm)
w=20;                                     % extent of mystacial pad M (mm)
Mh=0.5;                                   % Mass of hair (mg)
Mf=10;                                    % Mass of follicle (mg)
C=-1.43;                                  % Center of Mass (mm)
I0=112;                                   % Moment of inertia : (mg*(mm^2))

D=-3.5;                                     % Ld: length of muscle attachment from the skin (mm)
zeta_up=2.2;                                % upper dampers (mg/msec)
zeta_dn=2.2;                                % lower dampers (mg/msec)
zeta_gr=6;                               % vertical dampers (mg/msec)
Kup=0.2;                                  % upper springs constants (mg/msec) 
Kdn=0.2;                                  % Lower springs constants (mg/msec) 
Kgr=0.6;                                    % vertical springs constants (mg/msec)

%%% F-L curve
flow=0.67; 
fhigh=0.90;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Simulation time parameters (sec)
t_start=0;
t_step=0.001;
t_end=0.1;
time_shift=2*t_step;                        % Shifting left the output traces by xxx sec (time_shift=t_step, no shift)
%%%%%%%%%%%%%%%% END COMMON %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%   START LARGE ANGLES %%%%%%%%%%%%%%%
% %%% Large Angles
%     %%% Muscle Parameters 
r0 = 2.55;
tauc =7.4;  %msec
taur =5.0;  %msec
    
 %%%%%%%%%%%%%%%%%%% END  LARGE ANGLES %%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%  START UNITS TRANSLATION %%%%%%%%%%%%%%
% Match Units to MKS system 
Lf=Lf*1e-3;                                  
s=s*1e-3;                                    
w=w*1e-3; 

Mh=Mh*1e-6;                                
Mf=Mf*1e-6;                                  
C=C*1e-3; 
I0=I0*1e-12;


 tauc = tauc/1000;	
 taur =taur/1000;   
 %%% Bio-Mechanical Parametsrs (Pad)
 D=D/1000;     
 zeta_up=zeta_up/1000;  
 zeta_dn=zeta_dn/1000;
 zeta_gr=zeta_gr/1000;
 %%%%%%%%%%%%%%%%%%%  END UNITS TRANSLATION %%%%%%%%%%%%%%
 

