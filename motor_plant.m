% (c) Written By Erez Simony 2010, code for the model described in:  
% Simony, E., Bagdasarian K, Herfst L., Brecht M., Ahissar E, Golomb D. 
% Temporal and spatial characteristics of vibrissa responses to motor commands (2010). 
% Journal of Neuroscience, In press.

function [time_in_msec,delta_theta,delta_xc,delta_yc]=motor_plant( resting_angles, intrinsic_muscle_set, MN_spikes_times, force_factor)

% 
% General: the function simulate the spatio-temporal transformation between MN spikes into
% whiskers  movements, the transformation has 2 stages: from spikes to
% muscle force (CA2+ dynamics), and from muscle force to whisker movement (bio-mechanical model of the whisker pad)  
% The Bio-mechanical of the whisker pad is based on Hill et Al model (2008).

% The motor_plant function gets as an input <N>, the number of whiskers 
% connected in a single row, by (N+1) intrinsic muscles.
% <resting_angles> is a vector which define the resting angles of the whisker
% in a row (degree).
% Each intrinsic muscles can be activate or not (On, Off : '1','0')  by the
%  <intrinsic_muscle_set>  column vector. The first value ( muscle_set(1) ), refers to the most
% posterior intrinsic muscle. 
% <MN_spikes_times> is a time matrix , in which each row relates to a spike
% time series  (msec) from a single motoneuron (MN), which connects to a single intrinsic muscle. 
% The first row correspond to the most posterior intrinsic muscle .
% <force_factor> is the muscle force factor after the transformation, or
% the motor unit size . default is 1.


% global parameters
global N I C Lf M I0 Number_of_equations thetaRest K_springs  zeta_springs  rest_length 
global Xskin Yskin Xint Yint Xplate Yplate Xground Yground
global Xskin_dot Yskin_dot Xint_dot Yint_dot Xplate_dot Yplate_dot
global t_spikes Ci spk_idx  factor   
global spike_num f_enable time_shift  
global Xskin_anchors_right  Xskin_anchors_left Yskin_anchors Xplate_anchors_left Xplate_anchors_right Yplate_anchors D;
global zeta_up zeta_dn zeta_gr Kup Kdn Kgr 
global s w Mh Mf 
global t_start t_step t_end
global muscle_length_cont  muscle_length_cont_1 muscle_rest_length
%%%%%%%%%%%%%%%%      Bio-Mechnics of the whisker pad  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
Number_of_equations=6;
isi_train_1=1;
factor=force_factor; 

thetaRest=resting_angles;% 
I=[I0*ones(1,N) ];              
M=Mh+Mf;
equal_space=1;    % all springs are equal in length
space_factor=2;   % all springs are equal in length


t_length=length(t_start:t_step:t_end);
muscle_length_cont=[];
muscle_length_cont_1=[];



%%%%%%%%%%%%    Muscles Parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_enable=intrinsic_muscle_set;
if (factor==(length(isi_train_1)+1))
    t_spikes=MN_spikes_times*0.001;
else
    t_spikes=MN_spikes_times*0.001;
end

[int_num spike_num]=size(t_spikes);
Ci=zeros(N+1,spike_num);
spk_idx=f_enable;
%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% %%%


theta_vibrissa_init=[deg2rad(thetaRest) ]';
theta_vibrissa_vel_init=zeros(1,N);    
x_y_vibrissa_init=zeros(N,2);

        
if N>1
    vib_space_x=cumsum ( [(w-(N-1)*s)/2, s*ones(1,N-1) ]) ;
else
    vib_space_x=(w-(N-1)*s)/2;
end

vib_space_y=zeros(length(vib_space_x),1) ;
x_y_vibrissa_init=[ vib_space_x'  vib_space_y ];
x_vibrissa_vel_init=zeros(1,length(vib_space_x));         
y_vibrissa_vel_init=zeros(1,length(vib_space_y));        

x_pos_vel=[ vib_space_x; x_vibrissa_vel_init];
y_pos_vel=[vib_space_y';y_vibrissa_vel_init];
theta_pos_vel=[theta_vibrissa_init' ;theta_vibrissa_vel_init]; 
x_y_theta_pos_vel=[x_pos_vel;y_pos_vel;theta_pos_vel];


X1rest=x_y_vibrissa_init(1,1);
XNrest=x_y_vibrissa_init(N,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5  ANCORE POINTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % Constants boundaries  -  Posterior 
if equal_space==1
    Xskin_0fix=x_y_vibrissa_init(1,1)+C*cos(theta_vibrissa_init(1))-s/space_factor;
    Xplate_0fix=x_y_vibrissa_init(1,1)+(Lf+C)*cos(theta_vibrissa_init(1))-s/space_factor;
else
    Xskin_0fix=x_y_vibrissa_init(1,1)+C*cos(deg2rad(thetaRest(1)));
    Xplate_0fix=x_y_vibrissa_init(1,1)+(Lf+C)*cos(deg2rad(thetaRest(1)));
end
Yskin_0fix=-C*sin(deg2rad(thetaRest(1)));
Yplate_0fix=-(Lf+C)*sin(deg2rad(thetaRest(1)));
Xint_0fix=X1rest+C*cos(deg2rad(thetaRest(1)))-s;%X1rest+C*cos(deg2rad(thetaRest(1)))-2*s
Yint_0fix=-C*sin(deg2rad(thetaRest(1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants boundaries  -  Anterior 
if equal_space==1
    Xskin_sofix=x_y_vibrissa_init(N,1)+C*cos(theta_vibrissa_init(end))+s/space_factor;
    Xplate_sofix=x_y_vibrissa_init(N,1)+(Lf+C)*cos(theta_vibrissa_init(end))+s/space_factor;%w;
else
    Xskin_sofix=(x_y_vibrissa_init(N,1)+C*cos(deg2rad(thetaRest(end)))+w)/2;
    Xplate_sofix=(x_y_vibrissa_init(N,1)+(Lf+C)*cos(deg2rad(thetaRest(end)))+w)/2;%w;
end
Yskin_sofix=-C*sin(deg2rad(thetaRest(end)));
Yplate_sofix=-(Lf+C)*sin(deg2rad(thetaRest(end)));
Xint_sofix=XNrest+(Lf+C)*cos(deg2rad(thetaRest(end)))+2*s;%XNrest+(Lf+C)*cos(deg2rad(thetaRest(end)))+2*s;
Yint_sofix=-(Lf+C)*sin(deg2rad(thetaRest(end)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants boundaries  -  Vertical Springs 
Xground=[Xplate_0fix;x_y_vibrissa_init(:,1)+(Lf+C)*cos(deg2rad(thetaRest'));Xplate_sofix];
Yground=[Yplate_0fix;x_y_vibrissa_init(:,2)-(Lf+C)*sin(deg2rad(thetaRest'))-s/2;Yplate_sofix];


%%%%%%%%%%%%%%%%%%%%%%%%%%  R E ST POSITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% rest length of springs (first row- upper spring, second row -bottom spring)
rest_length=zeros(5,N);
Xskin=[Xskin_0fix;x_y_vibrissa_init(:,1)+C*cos(deg2rad(thetaRest'));Xskin_sofix];
Yskin=[Yskin_0fix;x_y_vibrissa_init(:,2)-C*sin(deg2rad(thetaRest'));Yskin_sofix];
Xplate=[Xplate_0fix;x_y_vibrissa_init(:,1)+(Lf+C)*cos(deg2rad(thetaRest'));Xplate_sofix];
Yplate=[Yplate_0fix;x_y_vibrissa_init(:,2)-(Lf+C)*sin(deg2rad(thetaRest'));Yplate_sofix];

%%%%%%%%%%%%%      R E S T     L E N G T H S                                            %%%%%%%%%%%%%%%%%%%%%
len_points=length(Xskin);

Xskin_anchors_left= Xskin(2:(len_points-1)) -s/space_factor;
Xskin_anchors_right=Xskin(2:(len_points-1)) +s/space_factor;
Yskin_anchors=Yskin(2:(len_points-1));

Xplate_anchors_left=Xplate(2:(len_points-1))-s/space_factor;
Xplate_anchors_right=Xplate(2:(len_points-1))+s/space_factor;
Yplate_anchors=Yplate(2:(len_points-1));



rest_length(1,:)= Xskin(2:(N+1),1)-Xskin_anchors_left;        % upper left
rest_length(2,:)= Xskin_anchors_right-Xskin(2:(N+1),1);       % upper right
rest_length(3,:)= Xplate(2:(N+1),1)-Xplate_anchors_left;      % bottom left
rest_length(4,:)= Xplate_anchors_right-Xplate(2:(N+1),1) ;    % bottom right
rest_length_ground=[sqrt( (Xplate-Xground).^2+(Yplate-Yground).^2 )]';
rest_length(5,:)=rest_length_ground(2:(end-1));

% SPRINGS CONSTANTS 
K_springs=zeros(5,N);
K_springs(1,:)=Kup*ones(1,N);                         % Upper Left Springs
K_springs(2,:)=Kup*ones(1,N);                         % Upper Right Springs
K_springs(3,:)=Kdn*ones(1,N);                         % Lower Left Springs
K_springs(4,:)=Kdn*ones(1,N);                         % Lower Right Springs
K_springs(5,:)=Kgr*ones(1,N);                         % Y-axis spring

%DAMPERS CONSTANTS
zeta_springs=zeros(5,N);
zeta_springs(1,:)=zeta_up*ones(1,N);                  % Upper Left dampers
zeta_springs(2,:)=zeta_up*ones(1,N);                  % Upper Right dampers
zeta_springs(3,:)=zeta_dn*ones(1,N);                  % Lower Left dampers
zeta_springs(4,:)=zeta_dn*ones(1,N);                  % Lower Right dampers
zeta_springs(5,:)=zeta_gr*ones(1,N);                  % Y-axis damper




%%%%%%%%%%%%%%%%%%%%%%%%%%       P O S I T I O N          %%%%%%%%%%%%%%%%%%%%%%%%%%%5
Xskin=[Xskin_0fix;x_y_vibrissa_init(:,1)+C*cos(theta_vibrissa_init);Xskin_sofix];
Yskin=[Yskin_0fix;x_y_vibrissa_init(:,2)-C*sin(theta_vibrissa_init);Yskin_sofix];
Xplate=[Xplate_0fix;x_y_vibrissa_init(:,1)+(Lf+C)*cos(theta_vibrissa_init);Xplate_sofix];
Yplate=[Yplate_0fix;x_y_vibrissa_init(:,2)-(Lf+C)*sin(theta_vibrissa_init);Yplate_sofix];
Xint=[Xint_0fix;x_y_vibrissa_init(:,1)-(D-C)*cos(theta_vibrissa_init);Xint_sofix];
Yint=[Yint_0fix;x_y_vibrissa_init(:,2)+(D-C)*sin(theta_vibrissa_init);Yint_sofix];
%%%%%%%%%%%%%%%%%%%%%%%%%%       V E L O C I T Y         %%%%%%%%%%%%%%%%%%%%%%%%%%%5
Xskin_dot=[0;x_vibrissa_vel_init'-C* theta_vibrissa_vel_init'.*sin(theta_vibrissa_init);0];
Yskin_dot=[0;-C* theta_vibrissa_vel_init'.*cos(theta_vibrissa_init);0];
Xplate_dot=[0;x_vibrissa_vel_init'-(Lf+C)* theta_vibrissa_vel_init'.*sin(theta_vibrissa_init);0];
Yplate_dot=[0;-(Lf+C)* theta_vibrissa_vel_init'.*sin(theta_vibrissa_init);0];
Xint_dot=[0;x_vibrissa_vel_init'+(D-C)*theta_vibrissa_vel_init'.*sin(theta_vibrissa_init);0];
Yint_dot=[0;(D-C)*theta_vibrissa_vel_init'.*cos(theta_vibrissa_init);0];

muscle_rest_length=[ sqrt( (Xint(2,1)-Xint(1,1) ).^2+ (Yint(2,1)-Yint(1,1) ).^2 ) ; sqrt( (Xint(3:(N+1),1)-Xskin(2:(N),1) ).^2+ (Yint(3:(N+1),1)-Yskin(2:(N),1) ).^2 ); sqrt( (Xint(end,1)-Xskin(end-1,1) ).^2+ (Yint(end,1)-Yskin(end-1,1) ).^2 ) ];
new_spk=t_spikes(intrinsic_muscle_set==1,:);
min_new_spk=min(new_spk(:,1));
tmp=new_spk==min_new_spk;
t_loop=[new_spk(tmp(1),:)];
t_loop_vir=[t_loop(1,:) (t_end+t_step)];

options = odeset('RelTol',[1e-6 ],'AbsTol',[1e-12 ]*ones(1,length(x_y_theta_pos_vel(:))));
x_y_init=x_y_theta_pos_vel(:);
T_time=zeros(t_length,1);
Y_SUM=zeros(t_length,N*6);
start=1;
for sim=1:length(t_loop)
    t_start=t_loop_vir(sim);
    t_end=t_loop_vir(sim+1);
    [T,Y] = ode45(@motor_plant_ode,[t_start:t_step:t_end],x_y_init,options);
    spk_idx= spk_idx+1;
    last=length(T)-1;
    T_time(start:(start+last-1))=T(1:last);     
    Y_SUM(start:(start+last-1),:)=Y(1:(last),:);
    start=start+last;
    x_y_init=[Y(end,:)]';
end
T=T_time;
Y=Y_SUM;


delta_xc=single(Y(:,1:Number_of_equations:Number_of_equations*N))- repmat(x_y_vibrissa_init(:,1)',length(T),1) ;
delta_xc=delta_xc((time_shift/t_step):end,:);
delta_yc=single(Y(:,3:Number_of_equations:Number_of_equations*N))- repmat(x_y_vibrissa_init(:,2)',length(T),1);    
delta_yc=delta_yc((time_shift/t_step):end,:);
delta_theta  =rad2deg(single(Y(:,5:Number_of_equations:Number_of_equations*N)))-repmat(thetaRest,length(T),1);
delta_theta=delta_theta((time_shift/t_step):end,:)-repmat(delta_theta((time_shift/t_step),:),length(T)-((time_shift/t_step)-1),1);
time_in_msec=T(1:(end-((time_shift/t_step)-1)))*1000;




