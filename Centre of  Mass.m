%BE2BM1
%Ekgari Kasawala

close all
clear

%% load files into matlab
StaticTF=importdata('4_Static_Trial_Forces.txt');                                    %static trial forces file
StaticTT=importdata('4_Static_Trial_Trajectories.txt');                              %static trial trajectory file
VerticalJF=importdata('4_Vertical_jump_Forces.txt');                                 %vertical jump forces file
VerticalJT=importdata('4_Vertical_jump_Trajectories.txt');  

%% Resampling  force files so that they have the same sampling rate as the trajectory files
% static Force file
OriginalSamplingF=1000;                                         %original sampling frequency
DesiredSamplingF=100;
[p,q]=rat(DesiredSamplingF/OriginalSamplingF);                            % this  finds a ration of the desired frequancy to the original one
RSTF=resample(StaticTF,p,q);                            %resampled static trial force file
% Vertival Jump force file
RVJF=resample(VerticalJF,p,q);

%% finding frame size and  total duration of each file
SF=100;
[f_maxS,~]=size(StaticTT);                                                      
Onestep=1/SF;                                    %time for one step
TotalTimeS=(1:f_maxS)*Onestep;                          %total duration for  static trial     
[f_maxV,~]=size(VerticalJT);   
TotalTimeV=(1:f_maxV)*Onestep;                           %total duration for vertcial jump trial
%% Static combined forces from static force file 
SF_combinedForces=RSTF(:,21:23);
%% Static trajectories for each segment
LASIS=StaticTT(:,3:5);
RASIS=StaticTT(:,6:8);
LPSIS=StaticTT(:,9:11);
RPSIS=StaticTT(:,12:14);
LEFT_THIGH=StaticTT(:,15:17);
LEFT_KNEE=StaticTT(:,18:20);
LEFT_TIBIA=StaticTT(:,21:23);
LEFT_ANKLE=StaticTT(:,24:26);
LEFT_HEEL=StaticTT(:,27:29);
LEFT_TOE=StaticTT(:,30:32);
RIGHT_THIGH=StaticTT(:,33:35);
RIGHT_KNEE=StaticTT(:,36:38);
RIGHT_TIBIA=StaticTT(:,39:41);
RIGHT_ANKLE=StaticTT(:,42:44);
RIGHT_HEEL=StaticTT(:,45:47);
RIGHT_TOE=StaticTT(:,48:50);

%% VERTICAL JUMP FORCES
VF1_Forces=RVJF(:,3:5);
VF1_Moment=RVJF(:,6:8);
VF1_COP=RVJF(:,9:11);
VF2_Force=RVJF(:,12:14);
VF2_Moment=RVJF(:,15:17);
VF2_COP=RVJF(:,18:20);
VF_combinedForces=RVJF(:,21:23);
VF_combinedMoments=RVJF(:,24:26);
VF_combinedCOP=RVJF(:,27:29);

%% VERTICAL JUMP TRAJECTORIES
VLASIS=VerticalJT(:,3:5);
VRASIS=VerticalJT(:,6:8);
VLPSIS=VerticalJT(:,9:11);
VRPSIS=VerticalJT(:,12:14);
VLEFT_THIGH=VerticalJT(:,15:17);
VLEFT_KNEE=VerticalJT(:,18:20);
VLEFT_TIBIA=VerticalJT(:,21:23);
VLEFT_ANKLE=VerticalJT(:,24:26);
VLEFT_HEEL=VerticalJT(:,27:29);
VLEFT_TOE=VerticalJT(:,30:32);
VRIGHT_THIGH=VerticalJT(:,33:35);
VRIGHT_KNEE=VerticalJT(:,36:38);
VRIGHT_TIBIA=VerticalJT(:,39:41);
VRIGHT_ANKLE=VerticalJT(:,42:44);
VRIGHT_HEEL=VerticalJT(:,45:47);
VRIGHT_TOE=VerticalJT(:,48:50);

%% constants   from Davis Hip model                                                           
LLegLength=sum(sqrt((LASIS(:,1)-LEFT_ANKLE(:,1)).^2+(LASIS(:,2)-LEFT_ANKLE(:,2)).^2+(LASIS(:,3)-LEFT_ANKLE(:,3)).^2))/289;                                                         %in mm
RLegLength=sum(sqrt((RASIS(:,1)-RIGHT_ANKLE(:,1)).^2+(RASIS(:,2)-RIGHT_ANKLE(:,2)).^2+(RASIS(:,3)-RIGHT_ANKLE(:,3)).^2))/289;                                                        %in mm
Weight=abs( sum(SF_combinedForces(:,3))/289);                                                               %in kgs
Rmarker=12;                                                                                           %in mm
SL=-1;     
SR=1;
LegLength=0.5*(LLegLength+RLegLength);
C=(0.115*LegLength)-0.0153;
Theta=(28.4*pi )/180;
Beta=(18*pi)/180;
Xdis=(0.1288*LegLength)-0.0153;
BAS=(0.1288*LegLength)-0.04856;
SubjectHeight=LegLength/0.530;    % from segment length chart
%% ASIS distance 
LA_RA=sqrt((VLASIS(:,1)-VRASIS(:,1)).^2+(VLASIS(:,2)-VRASIS(:,2)).^2+(VLASIS(:,3)-VRASIS(:,3)).^2);   % magnitude of Left asis minus right asis 
ASISDISTANCE=sum(LA_RA)/610;
OL= ((VLASIS-VRASIS) *0.5)+VRASIS;                               % origin of the local coordinate frame                                                                 %origin of the local coordinate frame 
GTL=OL';                                                          %translation from global to local frame
OLStatic=((LASIS-RASIS) *0.5)+RASIS; 

%% left and right hip coordinates in the local frame using Davis Hip Model

%Left Hip 
LHLX=((C*sin(Theta))-(0.5*ASISDISTANCE))*SL;
LHLY=((-Xdis-Rmarker)*cos(Beta))+(C*(cos(Theta)*sin(Beta)));
LHLZ=((-Xdis-Rmarker)*sin(Beta))-(C*(cos(Theta)*cos(Beta)));

%Right hip
RHLX=((C*sin(Theta))-(0.5*ASISDISTANCE))*SR;
RHLY=((-Xdis-Rmarker)*cos(Beta))+(C*(cos(Theta)*sin(Beta)));
RHLZ=((-Xdis-Rmarker)*sin(Beta))-(C*(cos(Theta)*cos(Beta)));

%% calculating GRL         
%first defining line 
I=(VRASIS-OL)./(sqrt((VRASIS(:,1)-OL(:,1)).^2+(VRASIS(:,2)-OL(:,2)).^2+(VRASIS(:,3)-OL(:,3)).^2)); %where i is vector in the x direction towards the Right ASIS 
i=I';

%2nd defining line 
L2=(VLPSIS-OL)./(sqrt((VLPSIS(:,1)-OL(:,1)).^2+(VLPSIS(:,2)-OL(:,2)).^2+(VLPSIS(:,3)-OL(:,3)).^2));     
K=(cross(I,L2));                                                                                                 
k=K';
J=(cross(K,I));                                                                                                  % where j is a vector component  in the y direction

j=J';
LH =[ LHLX LHLY LHLZ ]';     %location of left hip  in the local coordinate frame
RH=[RHLX RHLY RHLZ]';         %location of the right hip in the local coordinates 
    
for f=1:f_maxV
   GRL(1:3,1:3,f)=[i(1:3,f),j(1:3,f),k(1:3,f)];                            %rotation matrix                                                                  %global rotation to Local
   LHipG(1:3,f)=(GRL(:,:,f)*LH)+ GTL(:,f);                                 %location of left hip  in the global coordinate frame     
   RHipG(1:3,f)=(GRL(:,:,f)*RH)+ GTL(:,f);                                 %location of the right hip in the  global coordinate frame
end
% transpose of the left and right hip centres  back into column vectors
Hipright=RHipG';                     
Hipleft=LHipG';

for f=1:f_maxV         %% angle of excursions for left leg

    %segment angles using y(2) and z(3)
    %pelvic tilt
    PelvisAngleL(f)=rad2deg(atan2((VLPSIS(f,3)-VLASIS(f,3)),(VLPSIS(f,2)-VLASIS(f,2))))+90;    %left
    PelvisAngleR(f)=rad2deg(atan2((VRPSIS(f,3)-VRASIS(f,3)),(VRPSIS(f,2)-VRASIS(f,2))))+90;    %right
    %Thigh 
    ThighAngleL(f)=rad2deg(atan2((Hipleft(f,3)-VLEFT_KNEE(f,3)),(Hipleft(f,2)-VLEFT_KNEE(f,2))));
    ThighAngleR(f)=rad2deg(atan2((Hipright(f,3)-VRIGHT_KNEE(f,3)),(Hipright(f,2)-VRIGHT_KNEE(f,2))));
    %shank
    ShankAngleL(f)=rad2deg(atan2((VLEFT_KNEE(f,3)-VLEFT_ANKLE(f,3)),(VLEFT_KNEE(f,2)-VLEFT_ANKLE(f,2))));
    ShankAngleR(f)=rad2deg(atan2((VRIGHT_KNEE(f,3)-VRIGHT_ANKLE(f,3)),(VRIGHT_KNEE(f,2)-VRIGHT_ANKLE(f,2))));
    
    
    %foot
    FootAngleL(f)=rad2deg(atan2((abs(VLEFT_HEEL(f,3)-VLEFT_TOE(f,3))),(abs(VLEFT_HEEL(f,2)-(VLEFT_TOE(f,2))))));
    FootAngleR(f)=rad2deg(atan2((abs(VRIGHT_HEEL(f,3)-VRIGHT_TOE(f,3))),(abs(VRIGHT_HEEL(f,2)-(VRIGHT_TOE(f,2))))));
 % joint angles left leg
    KneeAngleL(f)=(ShankAngleL(f)-ThighAngleL(f));
    AnkleAngleL(f)=FootAngleL(f)-ShankAngleL(f)-90;
    HipAngleL(f)=(PelvisAngleL(f)-ThighAngleL(f));
  % joint angles right leg
    KneeAngleR(f)=(ShankAngleR(f)-ThighAngleR(f));
    AnkleAngleR(f)=FootAngleR(f)-ShankAngleR(f)-90;
    HipAngleR(f)=(PelvisAngleR(f)-ThighAngleR(f));
end

%% Plotting hip Angles
figure;
subplot(2,2,1)
plot(TotalTimeV,HipAngleL,'b')
hold on 
plot(TotalTimeV,HipAngleR,'r')
hold off 
legend('left hip','right hip','Location','Northwest')
title('HIP ANGLE OF EXCURSIONS')
xlabel('time(s)')
ylabel('Hip Angle (degrees)')

subplot(2,2,2)
plot(TotalTimeV,KneeAngleL,'b')
hold on 
plot(TotalTimeV,KneeAngleR,'r')
hold off
legend('left knee','right knee','Location','Northwest')
title('KNEE ANGLE OF EXCURSIONS')
xlabel('time(s)')
ylabel('Knee Angle (degrees)')

subplot(2,2,3)
plot(TotalTimeV,AnkleAngleL,'b')
hold on 
plot(TotalTimeV,AnkleAngleR,'r')
hold off
legend('left ankle','right ankle','Location','Northwest')
title('ANKLE ANGLE OF EXCURSIONS')
xlabel('time(s)')
ylabel('Ankle Angle (degrees)')
shg

%% Standard maximum height of centre of mass using the Pelvis origin
OLZ=OL(:,3);                                     % z coordinates of the pelvis origin
subplot(2,2,4)
plot(TotalTimeV,OLZ,'b')
hold on 
plot(TotalTimeS,OLStatic(:,3),'r')
hold off
title('Pelvis Origin Trajectory')
xlabel('y coordinate')
ylabel('Z coordinate')
HeightStandard=(max(OLZ)-max(OLStatic(:,3)))/1000;   %in m

%% Lower + Upper Body Segment COM method
%individual segment weight 
%using coeffients from dempster table tim;e the subject's body weight
ThighM=0.100*Weight;
ShankM=0.0465*Weight;
FootM=0.0145*Weight;
UpperBodyW=Weight-((ThighM+ShankM+FootM).*2);  %total body weight minus lower body weight

%COM positions each segment of  leftleg and right leg   in Vertical jump
%difference between joints on either end times the distal coefficients from the dempster table
LeftThighP=(0.567*(Hipleft-VLEFT_KNEE))+VLEFT_KNEE; 
LeftShankP=(0.567*(VLEFT_KNEE-VLEFT_ANKLE))+VLEFT_ANKLE;
LeftFootP=(0.50*(VLEFT_ANKLE-VLEFT_TOE))+VLEFT_TOE;
RightThighP=(0.567*(Hipright-VRIGHT_KNEE))+VRIGHT_KNEE; 
RightShankP=(0.567*(VRIGHT_KNEE-VRIGHT_ANKLE))+VRIGHT_ANKLE;
RightFootP=(0.50*(VRIGHT_ANKLE-VRIGHT_TOE))+VRIGHT_TOE;
UpperBodyP=OL;

%height calculation using the z component
%sum of com position of each segment times the weight divided by the sum of the weight
ZHeight=(((ThighM*LeftThighP(:,3))+(ShankM*LeftShankP(:,3))+(FootM*LeftFootP(:,3))+(ThighM*RightThighP(:,3))+(ShankM*RightShankP(:,3))+(FootM*RightFootP(:,3)+(UpperBodyP(:,3)*UpperBodyW))))/((2*(ThighM+ShankM+FootM))+UpperBodyW);
HeightMaxZ=(max(ZHeight)-ZHeight([1]))/1000;
figure ;
subplot (2,2,1)
plot(TotalTimeV,ZHeight)
title('ZCOM Height')
xlabel('time(s)')
ylabel('Height( m)')
%% kinematic  method
Mass=Weight/9.81;      % using W=mg
Force=VF_combinedForces(:,3); 
acceleration=(-Force-Weight)/Mass;  %F=ma

subplot (2,2,2)
plot(TotalTimeV,acceleration)
title('acceleration')
xlabel('time(s)')
ylabel('acceleration(ms.^-2)')
shg

velocity= cumtrapz(TotalTimeV,acceleration);   %integrates acceleration to find velocity
subplot (2,2,3)
plot(TotalTimeV,velocity,'r');
title('velocity')
xlabel('time(s)')
ylabel(' velocity(ms.^-1(')

displacement= cumtrapz(TotalTimeV,velocity); %integrates velocity to find displacement
subplot (2,2,4)
plot(TotalTimeV,displacement,'b');
title('displacement')
xlabel('time(s)')
ylabel('dispalcement(m)')

Heightkinematic=max(displacement);
%% projectile method
tO=find(max(velocity)==velocity)/100;  % intial time
th=find(min(velocity)==velocity)/100;  %final time
flighttime=(th-tO)/2 ;                 % duration of flight
g=9.81;
distance=((-g)*(flighttime.^2))/2;
HeightProjectile=abs(distance);

%% conservation of energy method
% a times t gives initial velocity since final velovity at maximum height is 0
initialvelocity=max(velocity);
HeightCOE=((initialvelocity).^2)/(2*g);

%% flight plot
figure ;
plot(TotalTimeV,acceleration,'b')
hold on
plot(TotalTimeV,velocity,'r')
hold off
hold on 
plot(TotalTimeV,displacement,'g')
hold off
xlabel('time(s)')
ylabel('Acceleration,velocity and displacement')
title('Flight Plot')
legend('acceleration','velocity','displacement','Location','Northwest')

%% Ground reaction plot
figure;
plot(TotalTimeV,-Force,'k')
xlabel('time(s)')
ylabel('Ground reaction')
title('GRF plot')
%% altanative calculation for Projectile method
%To get the first trough on the right
S=displacement;
maximumS=max(S);   %maximum displacement
Width=find(maximumS==S);
figure ;
plot(displacement,'b')
hold on 
plot(Width,maximumS,'k*')      % plots the maximum and corresponding  x value
hold off
index=Width;            %this is the start of the flight  wheree flight time is 0
widthr=0;           %distance from x value of the  peak to the right base 
while S(index+1)<S(index)
    widthr =widthr+1;
    index=index+1;
end
index=Width;% reset value of index back to peak  corresponding x value

%to get point on the left 
widthl=0;            % distance from x value of the peak to the base on the left
while S(index-1)<S(index)
    widthl=widthl+1;
    index=index-1;
end
hold on
plot(Width+widthr,S(Width+widthr),'ro')    % plots right trough
hold off
hold on
plot(Width-widthl,S(Width+widthl),'ro')
hold off
% base point plots are a bit off so will not be using this calculation
XT=Width-widthl:Width+widthr;
YV=S(Width-widthl:Width+widthr);
hold on
plot(XT,YV,'g')
hold   off

StartT=Width-widthl;
flightTime=((Width+widthr)-StartT)/100;
g=9.81;
PMHeight=g*(flightTime.^2)/8;   % too high 
