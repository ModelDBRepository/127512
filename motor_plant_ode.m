% (c) Written By Erez Simony 2010, code for the model described in:  
% Simony, E., Bagdasarian K, Herfst L., Brecht M., Ahissar E, Golomb D. 
% Temporal and spatial characteristics of vibrissa responses to motor commands (2010). 
% Journal of Neuroscience, In press.


function [dFdt] = motor_plant_ode(t,F)

 
global  N I C Lf M Number_of_equations  K_springs  zeta_springs  rest_length 
global Xskin Yskin Xint Yint Xplate Yplate Xground Yground
global Xskin_dot Yskin_dot Xint_dot Yint_dot Xplate_dot Yplate_dot
global r0 tauc taur muscle_rest_length flow fhigh
global t_spikes Ci spk_idx factor muscle_length_cont  muscle_length_cont_1
global f_enable 
global Xskin_anchors_left   Xskin_anchors_right Yskin_anchors Xplate_anchors_left Xplate_anchors_right Yplate_anchors D;
global  fl_stat

Number_of_equations=6;
dFdt = zeros(Number_of_equations*N,1);



Xc=F(1:Number_of_equations:Number_of_equations*N);
Xc_dot=F(2:Number_of_equations:Number_of_equations*N);
Yc=F(3:Number_of_equations:Number_of_equations*N);
Yc_dot=F(4:Number_of_equations:Number_of_equations*N);
theta=F(5:Number_of_equations:Number_of_equations*N);
theta_dot=F(6:Number_of_equations:Number_of_equations*N);


% attachment of position points  
Xskin(2:(N+1),1)=Xc+C*cos(theta);
Yskin(2:(N+1),1)=Yc-C*sin(theta);
Xint(2:(N+1),1)=Xc-(D-C)*cos(theta);
Yint(2:(N+1),1)=Yc+(D-C)*sin(theta);
Xplate(2:(N+1),1)=Xc+(Lf+C)*cos(theta);
Yplate(2:(N+1),1)=Yc-(Lf+C)*sin(theta);


% attachment velocity of position points  
Xskin_dot(2:(N+1),1)=Xc_dot-C*theta_dot.*sin(theta);
Yskin_dot(2:(N+1),1)=Yc_dot-C*theta_dot.*cos(theta);
Xint_dot(2:(N+1),1)=Xc_dot+(D-C)*theta_dot.*sin(theta);
Yint_dot(2:(N+1),1)=Yc_dot+(D-C)*theta_dot.*cos(theta);
Xplate_dot(2:(N+1),1)=Xc_dot-(Lf+C)*theta_dot.*sin(theta);
Yplate_dot(2:(N+1),1)=Yc_dot-(Lf+C)*theta_dot.*cos(theta);

% decoupled 
upper_left_lengths=sqrt(  (Xskin(2:(N+1),1)-Xskin_anchors_left ).^2+ (Yskin(2:(N+1),1)-Yskin_anchors ).^2);
upper_right_lengths=sqrt(  (Xskin(2:(N+1),1)-Xskin_anchors_right ).^2+ (Yskin(2:(N+1),1)-Yskin_anchors ).^2);
lower_left_lengths=sqrt( (Xplate(2:(N+1),1)-Xplate_anchors_left).^2+ (Yplate(2:(N+1),1)-Yplate_anchors ).^2 );
lower_right_lengths=sqrt( (Xplate(2:(N+1),1)-Xplate_anchors_right).^2+ (Yplate(2:(N+1),1)-Yplate_anchors ).^2 );
muscle_length=[ sqrt( (Xint(2,1)-Xint(1,1) ).^2+ (Yint(2,1)-Yint(1,1) ).^2 ) ; sqrt( (Xint(3:(N+1),1)-Xskin(2:(N),1) ).^2+ (Yint(3:(N+1),1)-Yskin(2:(N),1) ).^2 ); sqrt( (Xint(end,1)-Xskin(end-1,1) ).^2+ (Yint(end,1)-Yskin(end-1,1) ).^2 ) ];
ground_length=sqrt( (Xplate-Xground).^2+(Yplate-Yground).^2 );


for idx=0:N-1
    if (~isempty(t_spikes))
    
    if ( f_enable(idx+1) == 0)
         F_intrinsic=0;
    else
        sp_idx=spk_idx(idx+1,1);
        t_spk=t_spikes(idx+1, sp_idx);
     if (t>=t_spk)
        basic_exp=exp(-(t-t_spk)/tauc)-exp(-(t-t_spk)/taur);
           if ( sp_idx>1) 
                delta_t=t_spikes( idx+1,sp_idx)-t_spikes( idx+1,sp_idx-1);                                
                Ci( idx+1,sp_idx)=Ci( idx+1,sp_idx-1)*exp(-delta_t/tauc)+(r0*tauc/(tauc-taur))*( exp(-(delta_t)/tauc)-exp(-(delta_t)/taur) );
            end
         CA2=(r0*tauc/(tauc-taur))*basic_exp+Ci( idx+1,sp_idx)*exp(-(t-t_spk)/tauc);             
         muscle_length_cont_1(idx+1)=muscle_length(idx+1)/muscle_rest_length(idx+1);
         if (fl_stat)
             if (muscle_length_cont_1(idx+1)>=flow && muscle_length_cont_1(idx+1)<fhigh)
                    factor_length=((fhigh-flow)^-1)*(muscle_length_cont_1(idx+1)-flow);
              elseif muscle_length_cont_1(idx+1)< flow
                    factor_length=0;
             else
                    factor_length=1;
             end
         else
             factor_length=1;
         end
            
           F_intrinsic=((CA2.^4)./(1+CA2.^4))*1*factor(idx+1)*1e-3*factor_length;
     else
         F_intrinsic=0;

     end
     end
         
         if ( f_enable(idx+2) == 0)
            F_intrinsic_1=0;
        else
            sp_idx=spk_idx(idx+2,1);
            t_spk=t_spikes(idx+2, sp_idx);
        if (t>=t_spk) 
             basic_exp=exp(-(t-t_spk)/tauc)-exp(-(t-t_spk)/taur);
                if ( sp_idx>1) 
                    delta_t=t_spikes( idx+2,sp_idx)-t_spikes( idx+2,sp_idx-1);                                
                    Ci( idx+2,sp_idx)=Ci( idx+2,sp_idx-1)*exp(-delta_t/tauc)+(r0*tauc/(tauc-taur))*( exp(-(delta_t)/tauc)-exp(-(delta_t)/taur) );                                
                end
            CA2=(r0*tauc/(tauc-taur))*basic_exp+Ci( idx+2,sp_idx)*exp(-(t-t_spk)/tauc);      

            muscle_length_cont(idx+2)=muscle_length(idx+2)/muscle_rest_length(idx+2);
            if (fl_stat)
                if (muscle_length_cont(idx+2)>=flow && muscle_length_cont(idx+2)<fhigh)
                    factor_length=((fhigh-flow)^-1)*(muscle_length_cont(idx+2)-flow);
                elseif muscle_length_cont(idx+2)< flow
                    factor_length=0;
                else
                    factor_length=1;
                end
            else
                factor_length=1;
            end
                F_intrinsic_1=((CA2.^4)./(1+CA2.^4))*1*factor(idx+2)*1e-3*factor_length;
       else
         F_intrinsic_1=0;
       end
         end
         
         else
           F_intrinsic_1=0;
           F_intrinsic=0;
        end

         dFdt(1+idx*Number_of_equations) = F(2+idx*Number_of_equations);

                  
        % Lengths ....%%%%%%%%%%%%
        left_upper_length=upper_left_lengths(idx+1);
         right_upper_length=upper_right_lengths(idx+1);
        left_bottom_length=lower_left_lengths(idx+1);
        right_bottom_length=lower_right_lengths(idx+1);
        bottom_ground_length=ground_length(idx+2);
        
        % Velocity of "lengths"....%%%%%%%%%%%%
        left_upper_length_dot=(Xskin(idx+2)-Xskin_anchors_left(idx+1))*(Xskin_dot(idx+2)) +(Yskin(idx+2)-Yskin_anchors(idx+1))*(Yskin_dot(idx+2)) ;
        right_upper_length_dot=(Xskin_anchors_right(idx+1)-Xskin(idx+2))*(-Xskin_dot(idx+2)) +(Yskin_anchors(idx+1)-Yskin(idx+2))*(-Yskin_dot(idx+2)) ;
        left_bottom_length_dot=(Xplate(idx+2)-Xplate_anchors_left(idx+1))*(Xplate_dot(idx+2)) +(Yplate(idx+2)-Yplate_anchors(idx+1))*(Yplate_dot(idx+2)) ;
        right_bottom_length_dot=(Xplate_anchors_right(idx+1)-Xplate(idx+2))*(-Xplate_dot(idx+2)) +(Yplate_anchors(idx+1)-Yplate(idx+2))*(-Yplate_dot(idx+2)) ;

        bottom_ground_length_dot=(Xplate(idx+2)-Xground(idx+2))*(Xplate_dot(idx+2)) +(Yplate(idx+2)-Yground(idx+2))*(Yplate_dot(idx+2)) ;
        
         % Projection ....%%%%%%%%%%%%
        left_upper_proj_x=(Xskin(idx+2)-Xskin_anchors_left(idx+1))/left_upper_length;
         right_upper_proj_x=(Xskin_anchors_right(idx+1)-Xskin(idx+2))/right_upper_length;
        left_bottom_proj_x=(Xplate(idx+2)-Xplate_anchors_left(idx+1))/left_bottom_length;
        right_bottom_proj_x=(Xplate_anchors_right(idx+1)-Xplate(idx+2))/right_bottom_length;
        bottom_ground_proj_x=(Xplate(idx+2)-Xground(idx+2))/bottom_ground_length;
        
        left_upper_proj_y=(Yskin(idx+2)-Yskin_anchors(idx+1))/left_upper_length;
        right_upper_proj_y=(Yskin_anchors(idx+1)-Yskin(idx+2))/right_upper_length;
        left_bottom_proj_y=(Yplate(idx+2)-Yplate_anchors(idx+1))/left_bottom_length;
        right_bottom_proj_y=(Yplate_anchors(idx+1)-Yplate(idx+2))/right_bottom_length;
        bottom_ground_proj_y=(Yplate(idx+2)-Yground(idx+2))/bottom_ground_length;
                  
        if (idx==0)
            dFdt(2+idx*6)=(1/M)*(   -K_springs(1,1+idx)*( left_upper_length- rest_length(1,idx+1))*left_upper_proj_x+K_springs(2,1+idx)*( right_upper_length- rest_length(2,idx+1))*right_upper_proj_x...
            -K_springs(3,1+idx)*( left_bottom_length- rest_length(3,idx+1))*left_bottom_proj_x+K_springs(4,1+idx)*( right_bottom_length- rest_length(4,idx+1))*right_bottom_proj_x...
             -K_springs(5,1+idx)*(  bottom_ground_length- rest_length(5,idx+1))* bottom_ground_proj_x...
            +(-zeta_springs(1,1+idx)*( left_upper_length_dot/left_upper_length)*left_upper_proj_x+zeta_springs(2,1+idx)*( right_upper_length_dot/ right_upper_length)*right_upper_proj_x...
           -zeta_springs(3,1+idx)*( left_bottom_length_dot/left_bottom_length)*left_bottom_proj_x+zeta_springs(4,1+idx)*(  right_bottom_length_dot/ right_bottom_length)*right_bottom_proj_x)...
            -zeta_springs(5,1+idx)*( bottom_ground_length_dot/bottom_ground_length)*bottom_ground_proj_x...
            +(-F_intrinsic)*(Xint(idx+2)-Xint(idx+1))/sqrt( (Xint(idx+2)-Xint(idx+1))^2+(Yint(idx+2)-Yint(idx+1))^2 )...
           +(F_intrinsic_1)*(Xint(idx+3)-Xskin(idx+2))/sqrt( (Xint(idx+3)-Xskin(idx+2))^2+(Yint(idx+3)-Yskin(idx+2))^2 )    );
       
           
       else
                dFdt(2+idx*6)=(1/M)*(   -K_springs(1,1+idx)*( left_upper_length- rest_length(1,idx+1))*left_upper_proj_x+K_springs(2,1+idx)*( right_upper_length- rest_length(2,idx+1))*right_upper_proj_x...
            -K_springs(3,1+idx)*( left_bottom_length- rest_length(3,idx+1))*left_bottom_proj_x+K_springs(4,1+idx)*( right_bottom_length- rest_length(4,idx+1))*right_bottom_proj_x...
              -K_springs(5,1+idx)*(  bottom_ground_length- rest_length(5,idx+1))* bottom_ground_proj_x...
              +(-zeta_springs(1,1+idx)*( left_upper_length_dot/left_upper_length)*left_upper_proj_x+zeta_springs(2,1+idx)*( right_upper_length_dot/ right_upper_length)*right_upper_proj_x...
           -zeta_springs(3,1+idx)*( left_bottom_length_dot/left_bottom_length)*left_bottom_proj_x+zeta_springs(4,1+idx)*(  right_bottom_length_dot/ right_bottom_length)*right_bottom_proj_x)...
           -zeta_springs(5,1+idx)*( bottom_ground_length_dot/bottom_ground_length)*bottom_ground_proj_x...
           +(-F_intrinsic)*(Xint(idx+2)-Xskin(idx+1))/sqrt( (Xint(idx+2)-Xskin(idx+1))^2+(Yint(idx+2)-Yskin(idx+1))^2 )...
           +(F_intrinsic_1)*(Xint(idx+3)-Xskin(idx+2))/sqrt( (Xint(idx+3)-Xskin(idx+2))^2+(Yint(idx+3)-Yskin(idx+2))^2 )   );

        end
       

        dFdt(3+idx*Number_of_equations)=F(4+idx*Number_of_equations);
        if (idx==0)
             dFdt(4+idx*Number_of_equations)=(1/M)*(  -K_springs(1,1+idx)*( left_upper_length- rest_length(1,idx+1))*left_upper_proj_y+K_springs(2,1+idx)*( right_upper_length- rest_length(2,idx+1))*right_upper_proj_y...
            -K_springs(3,1+idx)*( left_bottom_length- rest_length(3,idx+1))*left_bottom_proj_y+K_springs(4,1+idx)*( right_bottom_length- rest_length(4,idx+1))*right_bottom_proj_y...
              -K_springs(5,1+idx)*(  bottom_ground_length- rest_length(5,idx+1))* bottom_ground_proj_y...
               -zeta_springs(5,1+idx)*( bottom_ground_length_dot/bottom_ground_length)*bottom_ground_proj_y...
            +(-zeta_springs(1,1+idx)*( left_upper_length_dot/left_upper_length)*left_upper_proj_y+zeta_springs(2,1+idx)*( right_upper_length_dot/ right_upper_length)*right_upper_proj_y...
            -zeta_springs(3,1+idx)*( left_bottom_length_dot/left_bottom_length)*left_bottom_proj_y+zeta_springs(4,1+idx)*(  right_bottom_length_dot/ right_bottom_length)*right_bottom_proj_y)...
            +(-F_intrinsic)*(Yint(idx+2)-Yint(idx+1))/sqrt( (Xint(idx+2)-Xint(idx+1))^2+(Yint(idx+2)-Yint(idx+1))^2 )...
           +(F_intrinsic_1)*(Yint(idx+3)-Yskin(idx+2))/sqrt( (Xint(idx+3)-Xskin(idx+2))^2+(Yint(idx+3)-Yskin(idx+2))^2 )  );
        else
             dFdt(4+idx*Number_of_equations)=(1/M)*(  -K_springs(1,1+idx)*( left_upper_length- rest_length(1,idx+1))*left_upper_proj_y+K_springs(2,1+idx)*( right_upper_length- rest_length(2,idx+1))*right_upper_proj_y...
            -K_springs(3,1+idx)*( left_bottom_length- rest_length(3,idx+1))*left_bottom_proj_y+K_springs(4,1+idx)*( right_bottom_length- rest_length(4,idx+1))*right_bottom_proj_y...
             -K_springs(5,1+idx)*(  bottom_ground_length- rest_length(5,idx+1))* bottom_ground_proj_y...
               -zeta_springs(5,1+idx)*( bottom_ground_length_dot/bottom_ground_length)*bottom_ground_proj_y...
            +(-zeta_springs(1,1+idx)*( left_upper_length_dot/left_upper_length)*left_upper_proj_y+zeta_springs(2,1+idx)*( right_upper_length_dot/ right_upper_length)*right_upper_proj_y...
            -zeta_springs(3,1+idx)*( left_bottom_length_dot/left_bottom_length)*left_bottom_proj_y+zeta_springs(4,1+idx)*(  right_bottom_length_dot/ right_bottom_length)*right_bottom_proj_y)...
            +(-F_intrinsic)*(Yint(idx+2)-Yskin(idx+1))/sqrt( (Xint(idx+2)-Xskin(idx+1))^2+(Yint(idx+2)-Yskin(idx+1))^2 )...
           +(F_intrinsic_1)*(Yint(idx+3)-Yskin(idx+2))/sqrt( (Xint(idx+3)-Xskin(idx+2))^2+(Yint(idx+3)-Yskin(idx+2))^2 )   );
        end
     
      
         dFdt(5+idx*Number_of_equations)=F(6+idx*Number_of_equations);
         
        
         
         if (idx==0)
             dFdt(6+idx*Number_of_equations)=(-1/I(idx+1))*( ((C)*cos(F(5+idx*Number_of_equations)))* (-K_springs(1,1+idx)*( left_upper_length- rest_length(1,idx+1))*left_upper_proj_y+K_springs(2,1+idx)*( right_upper_length- rest_length(2,idx+1))*right_upper_proj_y)...
                +(C*sin(F(5+idx*Number_of_equations)))*( -K_springs(1,1+idx)*( left_upper_length- rest_length(1,idx+1))*left_upper_proj_x+K_springs(2,1+idx)*( right_upper_length- rest_length(2,idx+1))*right_upper_proj_x)...
                +((Lf+C)*cos(F(5+idx*Number_of_equations)))*(-K_springs(3,1+idx)*( left_bottom_length- rest_length(3,idx+1))*left_bottom_proj_y+K_springs(4,1+idx)*( right_bottom_length- rest_length(4,idx+1))*right_bottom_proj_y)...
                +((Lf+C)*sin(F(5+idx*Number_of_equations)))*(-K_springs(3,1+idx)*( left_bottom_length- rest_length(3,idx+1))*left_bottom_proj_x+K_springs(4,1+idx)*( right_bottom_length- rest_length(4,idx+1))*right_bottom_proj_x)...
                +((Lf+C)*sin(F(5+idx*Number_of_equations)))*( -K_springs(5,1+idx)*(  bottom_ground_length- rest_length(5,idx+1))* bottom_ground_proj_x)...
                +((Lf+C)*cos(F(5+idx*Number_of_equations)))*( -K_springs(5,1+idx)*(  bottom_ground_length- rest_length(5,idx+1))* bottom_ground_proj_y)...
                +(C*cos(F(5+idx*Number_of_equations)))* (-zeta_springs(1,1+idx)*( left_upper_length_dot/left_upper_length)*left_upper_proj_y+zeta_springs(2,1+idx)*( right_upper_length_dot/ right_upper_length)*right_upper_proj_y)...
                +(C*sin(F(5+idx*Number_of_equations)))*( -zeta_springs(1,1+idx)*( left_upper_length_dot/left_upper_length)*left_upper_proj_x+zeta_springs(2,1+idx)*( right_upper_length_dot/ right_upper_length)*right_upper_proj_x)...
                +((Lf+C)*cos(F(5+idx*Number_of_equations)))*(-zeta_springs(3,1+idx)*( left_bottom_length_dot/left_bottom_length)*left_bottom_proj_y+zeta_springs(4,1+idx)*(  right_bottom_length_dot/ right_bottom_length)*right_bottom_proj_y)...
                +((Lf+C)*sin(F(5+idx*Number_of_equations)))*(-zeta_springs(3,1+idx)*( left_bottom_length_dot/left_bottom_length)*left_bottom_proj_x+zeta_springs(4,1+idx)*(  right_bottom_length_dot/ right_bottom_length)*right_bottom_proj_x) ...
                +((Lf+C)*cos(F(5+idx*Number_of_equations)))*(-zeta_springs(5,1+idx)*(  bottom_ground_length_dot/bottom_ground_length)*bottom_ground_proj_y)...
                +((Lf+C)*sin(F(5+idx*Number_of_equations)))*(-zeta_springs(5,1+idx)*(  bottom_ground_length_dot/bottom_ground_length)*bottom_ground_proj_x)...
                + ( cos(F(5+idx*Number_of_equations))*(  (D-C)*(F_intrinsic)*(Yint(idx+2)-Yint(idx+1))/sqrt( (Xint(idx+2)-Xint(idx+1))^2+(Yint(idx+2)-Yint(idx+1))^2 )... 
                + (C)*(F_intrinsic_1)*(Yint(idx+3)-Yskin(idx+2))/sqrt( (Xint(idx+3)-Xskin(idx+2))^2+(Yint(idx+3)-Yskin(idx+2))^2 ))...
                +sin(F(5+idx*Number_of_equations))*(  (D-C)*(F_intrinsic)*(Xint(idx+2)-Xint(idx+1))/sqrt( (Xint(idx+2)-Xint(idx+1))^2+(Yint(idx+2)-Yint(idx+1))^2 )... 
                + (C)*(F_intrinsic_1)*(Xint(idx+3)-Xskin(idx+2))/sqrt( (Xint(idx+3)-Xskin(idx+2))^2+(Yint(idx+3)-Yskin(idx+2))^2 )  ) )    );
            
                      
         else
               dFdt(6+idx*Number_of_equations)=(-1/I(idx+1))*( (C*cos(F(5+idx*Number_of_equations)))* (-K_springs(1,1+idx)*( left_upper_length- rest_length(1,idx+1))*left_upper_proj_y+K_springs(2,1+idx)*( right_upper_length- rest_length(2,idx+1))*right_upper_proj_y)...
                +(C*sin(F(5+idx*Number_of_equations)))*( -K_springs(1,1+idx)*( left_upper_length- rest_length(1,idx+1))*left_upper_proj_x+K_springs(2,1+idx)*( right_upper_length- rest_length(2,idx+1))*right_upper_proj_x)...
                +((Lf+C)*cos(F(5+idx*Number_of_equations)))*(-K_springs(3,1+idx)*( left_bottom_length- rest_length(3,idx+1))*left_bottom_proj_y+K_springs(4,1+idx)*( right_bottom_length- rest_length(4,idx+1))*right_bottom_proj_y)...
                +((Lf+C)*sin(F(5+idx*Number_of_equations)))*(-K_springs(3,1+idx)*( left_bottom_length- rest_length(3,idx+1))*left_bottom_proj_x+K_springs(4,1+idx)*( right_bottom_length- rest_length(4,idx+1))*right_bottom_proj_x) ...
                +((Lf+C)*sin(F(5+idx*Number_of_equations)))*( -K_springs(5,1+idx)*(  bottom_ground_length- rest_length(5,idx+1))* bottom_ground_proj_x)...
                +((Lf+C)*cos(F(5+idx*Number_of_equations)))*( -K_springs(5,1+idx)*(  bottom_ground_length- rest_length(5,idx+1))* bottom_ground_proj_y)...
                +((Lf+C)*cos(F(5+idx*Number_of_equations)))*(-zeta_springs(5,1+idx)*(  bottom_ground_length_dot/bottom_ground_length)*bottom_ground_proj_y)...
                +((Lf+C)*sin(F(5+idx*Number_of_equations)))*(-zeta_springs(5,1+idx)*(  bottom_ground_length_dot/bottom_ground_length)*bottom_ground_proj_x)...
                +(C*cos(F(5+idx*Number_of_equations)))* (-zeta_springs(1,1+idx)*( left_upper_length_dot/left_upper_length)*left_upper_proj_y+zeta_springs(2,1+idx)*( right_upper_length_dot/ right_upper_length)*right_upper_proj_y)...
                +(C*sin(F(5+idx*Number_of_equations)))*( -zeta_springs(1,1+idx)*( left_upper_length_dot/left_upper_length)*left_upper_proj_x+zeta_springs(2,1+idx)*( right_upper_length_dot/ right_upper_length)*right_upper_proj_x)...
                +((Lf+C)*cos(F(5+idx*Number_of_equations)))*(-zeta_springs(3,1+idx)*( left_bottom_length_dot/left_bottom_length)*left_bottom_proj_y+zeta_springs(4,1+idx)*(  right_bottom_length_dot/ right_bottom_length)*right_bottom_proj_y)...
                +((Lf+C)*sin(F(5+idx*Number_of_equations)))*(-zeta_springs(3,1+idx)*( left_bottom_length_dot/left_bottom_length)*left_bottom_proj_x+zeta_springs(4,1+idx)*(  right_bottom_length_dot/ right_bottom_length)*right_bottom_proj_x) ...
                + ( cos(F(5+idx*Number_of_equations))*(  (D-C)*(F_intrinsic)*(Yint(idx+2)-Yskin(idx+1))/sqrt( (Xint(idx+2)-Xskin(idx+1))^2+(Yint(idx+2)-Yskin(idx+1))^2 )... 
                + (C)*(F_intrinsic_1)*(Yint(idx+3)-Yskin(idx+2))/sqrt( (Xint(idx+3)-Xskin(idx+2))^2+(Yint(idx+3)-Yskin(idx+2))^2 ))...
                +sin(F(5+idx*Number_of_equations))*(  (D-C)*(F_intrinsic)*(Xint(idx+2)-Xskin(idx+1))/sqrt( (Xint(idx+2)-Xskin(idx+1))^2+(Yint(idx+2)-Yskin(idx+1))^2 )... 
                + (C)*(F_intrinsic_1)*(Xint(idx+3)-Xskin(idx+2))/sqrt( (Xint(idx+3)-Xskin(idx+2))^2+(Yint(idx+3)-Yskin(idx+2))^2 ) ) ));        
         end
         
     
                

end
end
