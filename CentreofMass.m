%BE2BM1
%Ekgari Kasawala

close all
clear
% load file into matlab
Data=importdata('Normal_Gait_Lower_Body(4).txt');
[f_max,~]=size(Data);                                                      
SamplingF=100;
Onestep=1/SamplingF;                                                       %time for onr step
TotalTime=(1:f_max)*Onestep;                                               %total duration of cycle


%% constants
Height=1837;                                                               %in mm
LeftLegLength=990;                                                         %in mm
RightLegLength=985;                                                        %in mm
Weght=70.2 ;                                                               %in kgs
Rmarker=12;                                                                %in mm
S=-1;                                                                      %for left leg
LegLength=0.5*(LeftLegLength+RightLegLength);
C=(0.115*LegLength)-0.0153;
Theta=(28.4*pi )/180;
Beta=(18*pi)/180;
Xdis=(0.1288*LegLength)-0.0153;
BAS= (0.1288*LegLength)-0.04856;                                           %Bilateral anteroposterior distance


%% Trajectories  for each marker
LeftASIS=Data(:,1:3);
RightASIS=Data(:,4:6);
LeftPSIS=Data(:,7:9);
RightPSIS=Data(:,10:12);
LeftThigh=Data(:,13:15);
LeftKnee=Data(:,16:18);
LeftTibia=Data(:,19:21);
Leftankle=Data(:,22:24);
LeftHeel=Data(:,25:27);
LeftToe=Data(:,28:30);
                                                                                                                                                                                          

%% finding ASIS distance
LEFTA_RIGHTA= sqrt((LeftASIS(:,1)-RightASIS(:,1)).^2+(LeftASIS(:,2)-RightASIS(:,2)).^2+(LeftASIS(:,3)-RightASIS(:,3)).^2); % magnitude of Left asis minus right asis
  ASIS_DISTANCE=sum(LEFTA_RIGHTA)/125;                                                                                     %average  asis distance

%% Hip joint centres calculation

%in the Local reference frame using the Davis Model
HLX=((C*sin(Theta))-(0.5*ASIS_DISTANCE))*S;
HLY=((-Xdis-Rmarker)*cos(Beta))+(C*(cos(Theta)*sin(Beta)));
HLZ=((-Xdis-Rmarker)*sin(Beta))-(C*(cos(Theta)*cos(Beta)));

% in the Global refrence frame using translation
OL= ((LeftASIS-RightASIS) *0.5)+RightASIS;                                                                     %origin of the local coordinate frame 
GTL=OL;                                                                                                        %translation from global to local frame
%first defining line 
i=(RightASIS-OL)./(sqrt((RightASIS(:,1)-OL(:,1)).^2+(RightASIS(:,2)-OL(:,2)).^2+(RightASIS(:,3)-OL(:,3)).^2)); %where i is vector in the x direction towards the Right ASIS 


%2nd defining line 
L2=(LeftPSIS-OL)./(sqrt((LeftPSIS(:,1)-OL(:,1)).^2+(LeftPSIS(:,2)-OL(:,2)).^2+(LeftPSIS(:,3)-OL(:,3)).^2));     
k=cross(i,L2);                                                                                                 

j=cross(k,i);                                                                                                  % where j is a vector component  in the y direction
 

    
 GRL=[ i(1:125)', j(1:125)', k(1:125)'];                                                                       %global rotation to Local
 %{
 GRL  is meant to be equal to  [ i  j k ] but was  returning as 125 by 9  . was unable to find a way to  make it  3
 by 3  with 125 frames .therefore i called first colum of each ( i j and k)
 and used them as my x y and z
 %}
 
 LH = [ HLX HLY HLZ ];                                                                                         %locaction of hip  in the local coordinate frame                                                                                                 %location of the hip in the global reference frame
 HipM=(GRL.*LH)+GTL;                                                                                           %location of him in the global   coordinate frame
             
 
 %% calculating  angle of excursions
 
%segment angles using x and z

for f=1:f_max
    
    
    %segment angles using x (1) and z(3)
    %pelvic tilt
   PelvisAngle(f)=rad2deg(atan2((LeftPSIS(f,3)-LeftASIS(f,3)),(LeftPSIS(f,1)-LeftASIS(f,1))))+90;
    
    %Thigh 
    ThighAngle(f)=rad2deg(atan2((HipM(f,3)-LeftKnee(f,3)),(HipM(f,1)-LeftKnee(f,1))));
    
    %shank
    ShankAngle(f)=rad2deg(atan2((LeftKnee(f,3)-Leftankle(f,3)),(LeftKnee(f,1)-Leftankle(f,1))));
    
    %foot
    FootAngle(f)=rad2deg(atan2((abs(LeftHeel(f,3)-LeftToe(f,3))),(abs(LeftHeel(f,1)-LeftToe(f,1)))));
    
 % joint angles
    KneeAngle(f)=ShankAngle(f)-ThighAngle(f);
    AnkleAngle(f)=FootAngle(f)-ShankAngle(f)-90;
    HipAngle(f)=(PelvisAngle(f)-ThighAngle(f));
end



%% plot Angle of excursions

F1=figure
subplot(2,2,1)
plot(TotalTime,HipAngle,'b')
title('HIP ANGLE OF EXCURSIONS')
xlabel('time(s)')
ylabel('Hip Angle (degrees)')

subplot(2,2,2)
plot(TotalTime,KneeAngle,'b')
title('KNEE ANGLE OF EXCURSIONS')
xlabel('time(s)')
ylabel('Knee Angle (degrees)')

subplot(2,2,3)
plot(TotalTime,AnkleAngle,'b')
title('ANKLE ANGLE OF EXCURSIONS')
xlabel('time(s)')
ylabel('Ankle Angle (degrees)')
shg

%% hip segment trajectory
%{
since  LeftPSIS and LeftASIS  were used to find  the pelvix tilt finding their
center of mass will give the pelvis  origin which can be used to  plot the
hip segment trajectory
%}
Coeficient=0.105;                                                              % used proximal coeficient from table in lecture slides
PCOM=LeftASIS+ (0.105*(LeftPSIS-LeftASIS));                                    %center of mass in the local coordinate frame
Pelvis_Origin=(GRL.*PCOM)+GTL;                                                 %in the global coorinate frame
x=Pelvis_Origin(:,1);                                                          %x coordinates
y=Pelvis_Origin(:,3);                                                          %y coordinates

%plot of the hip segment trajectory
subplot(2,2,4)
plot(x,y,'b')
title('HIP SEGMENT TRAJECTORY')
xlabel(' x coordnates')
ylabel('Y coordinates')
shg