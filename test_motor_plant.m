% (c) Written By Erez Simony 2010, code for the model described in:  
% Simony, E., Bagdasarian K, Herfst L., Brecht M., Ahissar E, Golomb D. 
% Temporal and spatial characteristics of vibrissa responses to motor commands (2010). 
% Journal of Neuroscience, In press.


global vib_num  resting_angles intrinsic_muscle_set force_factor  MN_spikes_times 
motor_plant_parameters_small_angles
% motor_plant_parameters_large_angles

% Call the motor_plant function
% Inputs: resting_angles, intrinsic_muscle_set,MN_spikes_times, force_factor
% Ouput:  time_in_msec,delta_theta,delta_xc,delta_yc 

[time_in_msec,delta_theta,delta_xc,delta_yc]=motor_plant(resting_angles, intrinsic_muscle_set, MN_spikes_times,force_factor);




% Plot whisker angle theta(degs) for "vib_num" and "vib_num-1" whiskers.
% (vib_num=1) , most posterior whisker.

figure

plot(time_in_msec,delta_theta(:,vib_num),'g','LineWidth',3)
% hold on
% plot(time_in_msec,delta_theta(:,vib_num-1),'k','LineWidth',3)
set(gca,'Position',[0.1759 0.1576 0.7705 0.7674],...
    'LineWidth',2,...
    'FontSize',16);
xlabel('Time (ms)','FontWeight','bold','FontSize',22);
ylabel('\theta (degs)','FontWeight','bold','FontSize',22);


% Plot whisker's center of mass translations Xc,Yc for "vib_num" 
figure
subplot(2,1,1,'LineWidth',2,'FontSize',16)
plot(time_in_msec,1000*delta_xc(:,vib_num),'g','LineWidth',3)
% hold on
% plot(time_in_msec,1000*delta_xc(:,vib_num-1),'k','LineWidth',3)
ylabel('x (mm)','FontSize',22,'FontName','Arial');

subplot(2,1,2,'LineWidth',2,'FontSize',16)
plot(time_in_msec,1000*delta_yc(:,vib_num),'g','LineWidth',3)
% hold on
% plot(time_in_msec,1000*delta_yc(:,vib_num-1),'k','LineWidth',3)
xlabel('Time (ms)','FontWeight','bold','FontSize',22);
ylabel('y (mm)','FontSize',22,'FontName','Arial');

