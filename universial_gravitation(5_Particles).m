%% CODE DESCRIPTION$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% code created by Jamie M Johns 2018 (tested in matlab 2013b) and can found at;
%  https://github.com/JamieMJohns/N-Particle-Simulation-of-Newton-s-Universal-Law-of-Gravitation-Matlab-

% code is created for demonstration of algorithm for Numerical application of Newton's universal law of gravitation


% Another matlab file,"universial_gravitation(N_Particles).m" (also found at github link [above]) demonstrates a
% more vectorized version of the algorithm and allow for the user to define
% N number of particles (planets) used in the simulation

%Sections of code:
%1 - 5 particle simulation in Two Dimensions
%2 - 5 particle simulation in Three Dimensions

%This code applies planetary data for sun,Merc,venus,earth,mars;
%->initial position (x,y,z)
%->mass of each planet
%->initial velocity (x,y,z)

%source of planetary data;
%  http://au.mathworks.com/help/physmod/sm/ug/model-planet-orbit-due-to-gravity.html




%% Section 1 - Example of algorithm for Universal gravitation in Two Dimensions

close all
clear all
clc
 
% initial position of each particle (planet)$$$$$$$$$$$$$$$$$$$
%e.g - p(5,:)=[x,y,z] position of particle 5 (units: meters)
p=[5.585e+08,5.585e+08;... %Sun (particle 1) position x,y
5.1979e+10,7.6928e+09;... %Merc (particle 2) position x,y
-1.5041e+10,9.708e+10;... %venus (particle 3) position x,y
-1.1506e+09,-1.391e+11;... %earth (particle 4) position x,y
-4.8883e+10,-1.9686e+11]; %mars (particle 5) position x,y


%mass of each particle (planet)
%e.g- m(4)=mass of 4th planet (units: kg)
m=[1.99e+30,3.3e+23,4.87e+24,5.97e+24,6.42e+23];

% initial velocity of each particle (planet)$$$$$$$$$$$$$$$$$$$
%e.g - v(2,:)=[x,y,z] velocity of particle 2 (units: meters
v=[-1.4663,11.124;... %Sun (particle 1) velocity x,y
-15205,44189;... %Merc (particle 2) velocity x,y
-34770,-5593.3;... %venus (particle 3) velocity x,y
29288,-398.58;... %earth (particle 4) velocity x,y
24533,-2762.2]; %mars (particle 5) velocity x,y


dt =24*60*60; %delta time (units: s); % difference in time between each instance of new position
              %[here I have used 24*60*60 which is number of seconds in one day]
G=6.673*10^(-11); %Universal Gravitational constant (units: m^3/(kg*s^2))
T=1000; % number of instances (time steps) to calculate new position
time=0; %recording of time for simulation time(1)=0
        %e.g= at instant 2, time(2)=time(1)+dt
        %and total simulated time will be equal to T*dt
        
        
p=repmat(p,1,1,T); % convert p position to 3D matrix where
                   % where p(3,:,53)=[x y z] position of particle 3
                   % at time step 53 (simulated time of time(t) [t=53])

%anonymous functions$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
R=@(pa,pb) sqrt((pa(1)-pb(1))^2+(pa(2)-pb(2))^2)+eps; %shortest distance between particle a and b
                                                                      %output is a scalar greater or equal than
                                                                      %machine epsilon (to avoid any potential division
                                                                      %of zero in function FG() )
r=@(pa,pb) pb-pa; %direction vector from particle b to particle a
FG=@(pa,pb,ma,mb) G.*ma.*mb.*r(pa,pb)./(R(pa,pb).^3); %equation of Force universal gravitation
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


pcm=0.0;%parameter used to show percentage of calculation complete
tic ;% start timer for length of simulation
commandwindow %bring up command window
fprintf('\n Now starting calculations:\n')

for t=2:T; % for instances (or timesteps) of 2 through to T (t=1 is the already defined initial state)
    for j=1:size(p,1); % for each particle j
        F=[0,0]; %initialise net force on particle j to be zero, F(1)=Fx F(2)=Fy F(3)=Fz
        for k=1:size(p,1); % for each particle k
            if j~=k %if j and k are not equal (not calculating a particles force on itself)
                p1=p(j,:,t-1); %get position of particle j [x,y,z] at previous instant (t-1)
                p2=p(k,:,t-1); %get position of particle k [x,y,z] at previous instant (t-1)
                F=F+FG(p1,p2,m(j),m(k)); %add gravitation for that particle k exerts on particle j                                     
            end                
        end
        v(j,:)=v(j,:)+dt.*F./m(j); %calculate velocity at current instant (t) for particle j
        p(j,:,t)=p(j,:,t-1)+dt.*v(j,:); %calculate and record position at current instant (t) for particle j
    end
time(t)=time(t-1)+dt; %determine time at instant t
%show progress of calculations
if t/T >= pcm; 
fprintf('%.0f%% of calculations complete at %.2fseconds\n',pcm*100,toc); %show    
pcm=pcm+0.1;
end

end
fprintf('Calculations complete at %.2fseconds',toc);
% visualisation###########################################################################################
hold on
cl=rand(size(p,1),1,3); %choose random color to plot each particle
                        %rand() is number between 0 and 1
                        % to define colour need 1x3 array = [R,G,B] (each between 0 and 1 (intensity)
sz=20*(1+1.*(m-min(m))./(max(m)-min(m))); %arbitrary equation to scale size of each plotted particle
                                            %based on particle mass (takes some experimenting to find something suitable)
t=1;                                            
for j=1:size(p,1); %plot each particle [plot at initial position]
%a seperate plot object is created for each particle
%to allow for each particle to have own size and color
plt(1,j)=plot(p(j,1,t),p(j,t),'b.','markers',sz(j),'color',cl(j,:)); %plot each particle with different size ('markers') and colour  
end
xlabel('x') %label x axis
ylabel('y') %label y axis
zlabel('z') %label z axis
title(sprintf('Time:%.0fseconds',0),'fontsize',15) %add title showing time (for the simulation) [to zero decimal places (%.0f)]
max3=@(x) max(max(max(abs(x))));
bnd=[-max3(p(:,1,:)) max3(p(:,1,:)) -max3(p(:,2,:)) max3(p(:,2,:))];

for t=1:1:T %for each time step
    for j=1:size(p,1); % for each particle
        set(plt(1,j),'Xdata',p(j,1,t),'Ydata',p(j,2,t))  %update x,y,z position of jth particle                                                                  
    end
title(sprintf('Time:%.0f days',time(t)./dt),'fontsize',15) %update title with latest time [to zero decimal places (%.0f)]
axis(bnd)

drawnow %update figure (visually) [above updates wont happen (visually) until you declare "drawnow"]
%dead=gsg
pause(0.005) %pause for 0.005 seconds (so animation is not too fast, you may wish to change pause time) 
end

% Note: you may wish to use the axis() function to scale axis of figure plot
% else, the axis in figure may constantly strecth/shrink as particles go beyond
% axis dimensions in subsequent time steps
%#################################################################################################################################

%% Section 2 - Example of algorithm for Universal gravitation in Three dimensions
close all
clear all
clc
 
% initial position of each particle (planet)$$$$$$$$$$$$$$$$$$$
%e.g - p(5,:)=[x,y,z] position of particle 5 (units: meters
p=[5.585e+08,5.585e+08,5.585e+08;... %Sun (particle 1) position x,y,z
5.1979e+10,7.6928e+09,-1.2845e+09;... %Merc (particle 2) position x,y,z
-1.5041e+10,9.708e+10,4.4635e+10;... %venus (particle 3) position x,y,z
-1.1506e+09,-1.391e+11,-6.033e+10;... %earth (particle 4) position x,y,z
-4.8883e+10,-1.9686e+11,-8.8994e+10]; %mars (particle 5) position x,y,z

%mass of each particle (planet)
%e.g- m(4)=mass of 4th planet (units: kg)
m=[1.99e+30,3.3e+23,4.87e+24,5.97e+24,6.42e+23];

% initial velocity of each particle (planet)$$$$$$$$$$$$$$$$$$$
%e.g - v(2,:)=[x,y,z] velocity of particle 2 (units: meters)
v=[-1.4663,11.124,4.837;... %Sun (particle 1) velocity x,y,z
-15205,44189,25180;... %Merc (particle 2) velocity x,y,z
-34770,-5593.3,-316.9;... %venus (particle 3) velocity x,y,z
29288,-398.58,-172.59;... %earth (particle 4) velocity x,y,z
24533,-2762.2,-1929.5]; %mars (particle 5) velocity x,y,z



dt =24*60*60; %delta time (units: s); % difference in time between each instance of new position
              %[here I have used 24*60*60 which is number of seconds in one day]
G=6.673*10^(-11); %Universal Gravitational constant (units: m^3/(kg*s^2))
T=1000; % number of instances (time steps) to calculate new position
time=0; %recording of time for simulation time(1)=0
        %e.g= at instant 2, time(2)=time(1)+dt
        %and total simulated time will be equal to T*dt
p=repmat(p,1,1,T); % convert p position to 3D matrix where
                   % where p(3,:,53)=[x y z] position of particle 3
                   % at time step 53 (simulated time of time(t) [t=53])

%anonymous functions$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
R=@(pa,pb) sqrt((pa(1)-pb(1))^2+(pa(2)-pb(2))^2+(pa(3)-pb(3))^2)+eps; %shortest distance between particle a and b
                                                                      %output is a scalar greater or equal than
                                                                      %machine epsilon (to avoid any potential division
                                                                      %of zero in function FG() )
r=@(pa,pb) pb-pa; %direction vector from particle b to particle a
FG=@(pa,pb,ma,mb) G.*ma.*mb.*r(pa,pb)./(R(pa,pb).^3); %equation of Force universal gravitation
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pcm=0.0;%parameter used to show percentage of calculation complete
tic ;% start timer for length of simulation
commandwindow %bring up command window
fprintf('\n Now starting calculations:\n')
% CALCULATIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
for t=2:T; % for instances (or timesteps) of 2 through to T (t=1 is the already defined initial state)
    for j=1:size(p,1); % for each particle j
        F=[0,0,0]; %initialise net force on particle j to be zero, F(1)=Fx F(2)=Fy F(3)=Fz
        for k=1:size(p,1); % for each particle k
            if j~=k %if j and k are not equal (not calculating a particles force on itself)
                p1=p(j,:,t-1); %get position of particle j [x,y,z] at previous instant (t-1)
                p2=p(k,:,t-1); %get position of particle k [x,y,z] at previous instant (t-1)
                F=F+FG(p1,p2,m(j),m(k)); %add gravitation for that particle k exerts on particle j                                     
            end                
        end
        v(j,:)=v(j,:)+dt.*F./m(j); %calculate velocity at current instant (t) for particle j
        p(j,:,t)=p(j,:,t-1)+dt.*v(j,:); %calculate and record position at current instant (t) for particle j
    end
    time(t)=time(t-1)+dt; %determine time at instant t

%show progress of calculations
if t/T >= pcm; 
fprintf('%.0f%% of calculations complete at %.2fseconds\n',pcm*100,toc); %show    
pcm=pcm+0.1;
end

end
fprintf('Calculations complete at %.2fseconds',toc);
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


% visualisation###########################################################################################
hold on
cl=rand(size(p,1),1,3); %choose random color to plot each particle
                        %rand() is number between 0 and 1
                        % to define colour need 1x3 array = [R,G,B] (each between 0 and 1 (intensity)
sz=20*(1+1.*(m-min(m))./(max(m)-min(m))); %arbitrary equation to scale size of each plotted particle
                                            %based on particle mass (takes some experimenting to find something suitable)
t=1;  

%optional (for setting boundary of x,y,z axis in plotted figure)##############
max3=@(x) max(max(max(x)));
min3=@(x) min(min(min(x)));
bnd=[min3(p(:,1,:)) max3(p(:,1,:)) min3(p(:,2,:)) max3(p(:,2,:)) min3(p(:,3,:)) max3(p(:,3,:))]; % get max/min values of position (x,y,z)
%###########################################################################

for j=1:size(p,1); %plot each particle [plot at initial position]
%a seperate plot object is created for each particle
%to allow for each particle to have own size and color
plt(1,j)=plot3(p(j,1,t),p(j,t),p(j,3,t),'b.','markers',sz(j),'color',cl(j,:)); %plot each particle with different size ('markers') and colour  
end
xlabel('x') %label x axis
ylabel('y') %label y axis
zlabel('z') %label z axis
title(sprintf('Time:%.4fseconds',0),'fontsize',15) %add title showing time (for the simulation) [to zero decimal places (%.0f)]  


for t=1:1:T %for each time step
    for j=1:size(p,1); % for each particle
        set(plt(1,j),'Xdata',p(j,1,t),'Ydata',p(j,2,t),'Zdata',p(j,3,t))  %update title with latest time [to zero decimal places (%.0f)]                                                               
    end
title(sprintf('Time:%.4f days',time(t)./dt),'fontsize',15) %update title with latest time [to zero decimal places (%.0f)]
view(3) %set figure view to 3 dimensions
axis(bnd)
drawnow %update figure (visually) [above updates wont happen (visually) until you declare "drawnow"]
pause(0.005) %pause for 0.005 seconds (so animation is not too fast, you may wish to change pause time) 
end

% Note: you may wish to use the axis() function to scale axis of figure plot
% else, the axis in figure may constantly strecth/shrink as particles go beyond
% axis dimensions in subsequent time steps
%#################################################################################################################################



