%% CODE DESCRIPTION$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% code created by Jamie M Johns 2018 (tested in matlab 2013b) and can found at;
%  https://github.com/JamieMJohns/N-Particle-Simulation-of-Newton-s-Universal-Law-of-Gravitation-Matlab-

% code is created for demonstration of Numerical application of Newton's universal law of gravitation
% code is also subject to improvement and optimization.

%Sections of code:
%1 - N particle simulation in Two Dimensions
%2 - N particle simulation in Three Dimensions


%This code applies same data as universal_gravitation.m for inital position (x,y,z) and mass
%as a bounds for generating N random particles by position (x,y,z) and mass.
%(initial velocities are set to zero for all particles).

%Other features of the code ("universial_gravitation(5_Particles).m");
%  -> calculation of force and new velocity/position is partially vectorised (improving speed of plotting) 
%  -> particles are bounded by maximum x,y,z;
%     if particles go outside of bounds their position is readjusted and velocity reset to zero
%  -> each instance of new particle position is plotted during calculation
%  ->since particles are plotted after each instant, positions are not stored in 3D matrix
%    but are instead kept as 2D matrix and overwritten at each calculation


%source of planetary data;
%  http://au.mathworks.com/help/physmod/sm/ug/model-planet-orbit-due-to-gravity.html

%% Section 1 - N particle simulation in Two Dimensions


close all %close all figure
clear all %clear all variables
clc %clear command window
 
N=1000; %number of particles to simulat
dt =24*60*60; %delta time (units: s); % difference in time between each instance of new position
              %[here I have used 24*60*60 which is number of seconds in one day]
T=300; % number of instances (time steps) to calculate new position
G=6.673*10^(-11); %Universal Gravitational constant (units: m^3/(kg*s^2))
time=0; %recording of time for simulation time(1)=0
        %e.g= at instant 2, time(2)=time(1)+dt
        %and total simulated time will be equal to T*dt
        

%sample initial positions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%po (below) will be used to determine a maximum and minimum x,y,z
po=[5.585e+08,5.585e+08;... %Sun position x,y
5.1979e+10,7.6928e+09;... %Merc position x,y
-1.5041e+10,9.708e+10;... %venus position x,y
-1.1506e+09,-1.391e+11;... %earth  position x,y
-4.8883e+10,-1.9686e+11]; %mars position x,y

po=po*10; %multiplied positions for scaling max/min x,y,z

%get N random X,Y coordinates (1 for each particle);
%where;  min(po)< X < max(po)
%and;  min(po)< Y < max(po)
minxy=min(min(po)); %mininum x,y,z
maxxy=max(max(po)); %maximum x,y,z
p=[minxy+rand(N,2).*(maxxy-minxy)]; % generate random position of N particles between minxy and maxxy
bnd=[minxy maxxy minxy maxxy]; %define bounds of position = [min(x) max(x) min(y) max(y)];
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%Sample random mass!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mo=[1.99e+30,3.3e+23,4.87e+24,5.97e+24,6.42e+23]; %mass of sun,Merc,venus,earth,mars
%get N random mass (1 for each particle)
%where;   min(mo)<mass<max(mo)
m=min(mo)+rand(1,N).*(max(mo)-min(mo));
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% initial velocity of each particle (planet)
v=zeros(N,2); %let velocity x,y equal zero initially for all particles


        
%anonymous functions$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
R=@(pa,pb) repmat(sqrt((pa(:,1)-pb(:,1)).^2+(pa(:,2)-pb(:,2)).^2)+eps,1,2); %shortest distance between particle a and b
                                                                      %output is a scalar greater or equal than
                                                                      %machine epsilon (to avoid any potential division
                                                                      %of zero in function FG() )
r=@(pa,pb) pb-pa; %direction vector from particle b to particle a
FG=@(pa,pb,ma,mb) G.*ma.*mb.*r(pa,pb)./(R(pa,pb).^3); %equation of Force universal gravitation
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


% INITIAL PLOT OF PARTICLES ###################################################
figure %create figure
hold on %keep each plotted object without deletion
cl=rand(size(p,1),1,3); %choose random color to plot each particle
                        %rand() is number between 0 and 1
                        % to define colour need 1x3 array = [R,G,B] (each between 0 and 1 (intensity)
sz=20*(1+1.*(m-min(m))./(max(m)-min(m))); %arbitrary equation to scale size of each plotted particle
                                            %based on particle mass (takes some experimenting to find something suitable)                                          
for j=1:size(p,1); %plot each particle [plot at initial position]
%a seperate plot object is created for each particle
%to allow for each particle to have own size and color
plt(1,j)=plot(p(j,1),p(j,2),'b.','markers',sz(j),'color',cl(j,:)); %plot each particle with different size ('markers') and colour  
end
xlabel('x') %label x axis
ylabel('y') %label y axis
zlabel('z') %label z axis
title(sprintf('Time:%.4fseconds',0),'fontsize',15) %add title showing time (for the simulation)

%################################################################################



k=1:size(p,1); %list of indexes for each planet

for t=1:T; % for T time steps
    p0=p;%position at initial timestep t-1
    v0=v;%velocity at initial timestep t-1
    for j=1:size(p,1); % for each particle j
        pj=repmat(p0(j,:),size(p0,1)-1,1);  % get position of particle j
                                          % and represent as matrix of size (N-1)x2
                                          %[column 1=x column 2=y]
        pk=p0(k~=j,:); %get position of all particles that are not j
                      %size=(N-1)x2  [size(pk)=size(pj)]
        massj=repmat(m(j),size(pj)); %get mass j and represent as (N-1)x2 matrix [same size as pj]
        massk=repmat(m(k~=j)',1,2); %get all mass that not j  and represent as (N-1)x2 matrix [same size as pj]
        F=sum(FG(pj,pk,massj,massk),1); %calculate net force                                
        v(j,:)=v0(j,:)+dt.*F./m(j); %calculate velocity (x,y) at current instant (t) for particle j
        p(j,:)=p0(j,:)+dt.*v(j,:); %calculate new position (x,y) for particle j at current instant t
        %check if new particle j position in boundary$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        % if no, set velocity to zero (inelastic collision)
        v(j,:)=((p(j,1))>bnd(1))*((p(j,1))<bnd(2)).*((p(j,2))>bnd(3))*((p(j,2))<bnd(4)).*v(j,:); %if particle i
        % if no for position x, set position x to be with in bounds
        p(j,1)=((p(j,1))<bnd(1)).*bnd(1)+((p(j,1))>bnd(2)).*bnd(2)+((p(j,1))>bnd(1))*((p(j,1))<bnd(2)).*p(j,1);
        % if no for position y, set position y to be with in bounds
        p(j,2)=((p(j,2))<bnd(3)).*bnd(3)+((p(j,2))>bnd(4)).*bnd(4)+((p(j,2))>bnd(3))*((p(j,2))<bnd(4)).*p(j,2);
        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    end
time=time+dt; %determine time at instant t

% UPDATE FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:size(p,1); % for each particle
        set(plt(1,j),'Xdata',p(j,1),'Ydata',p(j,2))  %update x,y,z position of jth particle                                                                  
    end
title(sprintf('Time:%.0f days',time./dt),'fontsize',15) %update title with latest time [to zero decimal places (%.0f)]
axis(bnd)

drawnow %update figure (visually) [above updates wont happen (visually) until you declare "drawnow"]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



%% Section 2 - N particle simulation in Three Dimensions


close all %close all figure
clear all %clear all variables
clc %clear command window
 
N=100; %number of particles to simulat
dt =24*60*60/2; %delta time (units: s); % difference in time between each instance of new position
              %[here I have used 24*60*60 which is number of seconds in one day]
T=1000; % number of instances (time steps) to calculate new position
G=6.673*10^(-11); %Universal Gravitational constant (units: m^3/(kg*s^2))
time=0; %recording of time for simulation time(1)=0
        %e.g= at instant 2, time(2)=time(1)+dt
        %and total simulated time will be equal to T*dt
        

%sample initial positions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%po (below) will be used to determine a maximum and minimum x,y,z
po=[5.585e+08,5.585e+08,5.585e+08;... %Sun  position x,y,z
5.1979e+10,7.6928e+09,-1.2845e+09;... %Merc position x,y,z
-1.5041e+10,9.708e+10,4.4635e+10;... %venus position x,y,z
-1.1506e+09,-1.391e+11,-6.033e+10;... %earth position x,y,z
-4.8883e+10,-1.9686e+11,-8.8994e+10]; %mars position x,y,z


po=po.*2; %multiplied positions for scaling max/min x,y,z

%get N random X,Y coordinates (1 for each particle);
%where;  min(po)< X < max(po)
%and;  min(po)< Y < max(po)
minxyz=min(min(min(po))); %mininum x,y,z
maxxyz=max(max(max(po))); %maximum x,y,z
p=[minxyz+rand(N,3).*(maxxyz-minxyz)]; % generate random position of N particles between minxy and maxxy
bnd=[minxyz maxxyz minxyz maxxyz minxyz maxxyz]; %define bounds of position = [min(x) max(x) min(y) max(y) min(z) max(z)];
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%Sample random mass!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mo=[1.99e+30,3.3e+23,4.87e+24,5.97e+24,6.42e+23]; %mass of sun,Merc,venus,earth,mars
%get N random mass (1 for each particle)
%where;   min(mo)<mass<max(mo)
m=min(mo)+rand(1,N).*(max(mo)-min(mo));
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% initial velocity of each particle (planet)
v=zeros(N,3); %let velocity x,y ,z equal zero initially for all particles


        
%anonymous functions$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
R=@(pa,pb) repmat(sqrt((pa(:,1)-pb(:,1)).^2+(pa(:,2)-pb(:,2)).^2+(pa(:,3)-pb(:,3)).^2)+eps,1,3); %shortest distance between particle a and b
                                                                      %output is a scalar greater or equal than
                                                                      %machine epsilon (to avoid any potential division
                                                                      %of zero in function FG() )
r=@(pa,pb) pb-pa; %direction vector from particle b to particle a
FG=@(pa,pb,ma,mb) G.*ma.*mb.*r(pa,pb)./(R(pa,pb).^3); %equation of Force universal gravitation
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


% INITIAL PLOT OF PARTICLES ###################################################
figure %create figure
hold on %keep each plotted object without deletion
cl=rand(size(p,1),1,3); %choose random color to plot each particle
                        %rand() is number between 0 and 1
                        % to define colour need 1x3 array = [R,G,B] (each between 0 and 1 (intensity)
sz=10*(1+2.*(m-min(m))./(max(m)-min(m))); %arbitrary equation to scale size of each plotted particle
                                            %based on particle mass (takes some experimenting to find something suitable)                                          
for j=1:size(p,1); %plot each particle [plot at initial position]
%a seperate plot object is created for each particle
%to allow for each particle to have own size and color
plt(1,j)=plot3(p(j,1),p(j,2),p(j,3),'b.','markers',sz(j),'color',cl(j,:)); %plot each particle with different size ('markers') and colour  
end
xlabel('x') %label x axis
ylabel('y') %label y axis
zlabel('z') %label z axis
title(sprintf('Time:%.4fseconds',0),'fontsize',15) %add title showing time (for the simulation)

%################################################################################



k=1:size(p,1); %list of indexes for each planet
for t=1:T; % for T time steps
    p0=p;%position at initial timestep t-1
    v0=v;%velocity at initial timestep t-1    
    for j=1:size(p,1); % for each particle j
        pj=repmat(p0(j,:),size(p0,1)-1,1);  % get position of particle j
                                          % and represent as matrix of size (N-1)x2
                                          %[column 1=x column 2=y]
        pk=p0(k~=j,:); %get position of all particles that are not j
                      %size=(N-1)x2  [size(pk)=size(pj)]
        massj=repmat(m(j),size(pj)); %get mass j and represent as (N-1)x2 matrix [same size as pj]
        massk=repmat(m(k~=j)',1,3); %get all mass that not j  and represent as (N-1)x2 matrix [same size as pj]
        F=sum(FG(pj,pk,massj,massk),1); %calculate net force                                
        v(j,:)=v0(j,:)+dt.*F./m(j); %calculate velocity (x,y) at current instant (t) for particle j
        p(j,:)=p0(j,:)+dt.*v(j,:); %calculate new position (x,y) for particle j
        %check if new particle j position in boundary$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        % if no, set velocity to zero (inelastic collision)
        v(j,:)=((p(j,1))>bnd(1))*((p(j,1))<bnd(2)).*((p(j,2))>bnd(3))*((p(j,2))<bnd(4)).*((p(j,3))>bnd(5)).*((p(j,3))<bnd(6)).*v(j,:);
        % if no for position y, set position x to be with in bounds
        p(j,1)=((p(j,1))<bnd(1)).*bnd(1)+((p(j,1))>bnd(2)).*bnd(2)+((p(j,1))>bnd(1))*((p(j,1))<bnd(2)).*p(j,1);
        % if no for position y, set position y to be with in bounds
        p(j,2)=((p(j,2))<bnd(3)).*bnd(3)+((p(j,2))>bnd(4)).*bnd(4)+((p(j,2))>bnd(3))*((p(j,2))<bnd(4)).*p(j,2);
        % if no for position y, set position z to be with in bounds
        p(j,3)=((p(j,3))<bnd(5)).*bnd(5)+((p(j,3))>bnd(6)).*bnd(6)+((p(j,3))>bnd(5))*((p(j,3))<bnd(6)).*p(j,3);        
        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    end
time=time+dt; %determine time at instant t

% UPDATE FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:size(p,1); % for each particle
        set(plt(1,j),'Xdata',p(j,1),'Ydata',p(j,2),'Zdata',p(j,3))  %update x,y,z position of jth particle                                                                  
    end
title(sprintf('Time:%.0f days',time./dt),'fontsize',15) %update title with latest time [to zero decimal places (%.0f)]
axis(bnd)
view(3)
drawnow %update figure (visually) [above updates wont happen (visually) until you declare "drawnow"]
pause(0.001)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
