%% universal Gravity N-particle run codestart=0 in command first
%code created by Jamie M Johns 2015  [tested with Matlab version R2013b]

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%PRECAUTION: BEFORE RUNNING THIS CODE MAKE SURE THAT ALL OPEN FILES ON YOUR
%            COMPUTER ARE SAVED (REGARDLESS OF WHETHER THE CODE HAS BEEN
%            MODIFIED FROM DEFAULT DOWLOADED CODE).
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 %     note: This code was created for fun and is subject to improvement and/or optimisation.
 
 %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        close all
        clear all
        clc  
      tic  
%parameters to modify#########################################################################################################

% General parameters  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Maxmass=20.^8; % maximum "particle" mass (kg)[units:kg]
        Minmass=20.^5; % mininum "particle" mass (kg) [units:kg]
        partnum=10000; % total number of "particles"  (N)
                     % each particle is given a random initial position (x,y,z) and random mass (between Minmass and Maxmass (defined above))        

        G=6.674.*(10^-5); %gravitational constant (default(theoretical) is approx: G =6.674.*(10^-11) [units: N*((m^2)/(kg^2))] )
        domainwidth=10000;  %length of sides the square area/volume that the "particles" are bounded within [meters]       
        threedim=2; %1=three dimension motion/space , other=two dimension motion/space 
        rotrate3d=1;  % 3d rotations of domain per 1 minute of footage [applies to 3D space/motion rendering]
        
        elasticcoll=1; % 1=particles rebound elastically from domain walls, other=particles don't rebound from wall
        
        ECP=0; % for elasticcoll=1 ; ECP=proportion of incident velocity that 
                  % contributes to particle rebound velocity for particle that
                  % has elastic collision with wall. [ECP>0 , range between 0 and 1] 
                  % (i.e rebound velocity = ECP*-1*(incident velocity) [direction is reversed after collision hence -1]
                  
         showdetail=1; %if showdetail=1; details of set parameters will be displayed at beginning of rendered output (useful with making video recordings)
         
         detailtime=1; % time,in seconds, that details of simulation are displayed at beginning of output (only if showdetail=1);
                        %  min is 1 second[detailtime=1];
                        
         mindist=1;  %mindist sets a mininum distance that can be detected between two particles;
                    %overall-->stops matlab from registering inf or nan for forces
                    %          and velcoties of particles that have seperation distance approaching zero meters

                    % you could set mindist to a value depending on value of "domainwidth above"
                    % examples -  for  domainwidth=2000 , i used mindist=1  (but you can ofcourse try something like mindist=0.000000001)
                    %             for  domainwidth=10  , i may use mindist=0.00000000001

                    % Setting mindist to a value not equal to zero avoids the problem in the applied set of equations;            
                        % i.e-  if, d12=0 (absolute seperation distance between particle 1 and particle 2)
                        %  where, F=-G*m1*m2/(d12^2)    [F=universial gravity force]
                        %  then, F->infinity  as d12->0 
                        %  as, F=-G*m1*m2/(0)=infinity (or nan)
                        %  hence,  velocities->inifinity
                        %  and so,   if mindist=C   then d12=(d12<=C)*C+(d12>C)*d12
                        %  therefore, maximum force possible: F--->-G*m1*m2/(C^2) [and not infinity]



                  
        outputpossize=[ 860, 100, 1049, 895];  %outputpossize=0--->default resolution and position of output
        %outputpossize=[ Px, Py, Ox, Oy]; %[Px Py]= position of output with respect to bottom-left corner of monitor 
                                          %[Ox Oy]=[width height] of the output
                                          
        %Example:  outputpossize=[ 860, 100, 1049, 895];  <-----good for 1080p monitor
        fps=30;   % frames per second for render      
        
        timelength=90;  % desired timelength for render [total frames to render=timelength*fps]

 % video recording parameters  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        record=1; % if record=1--->record output as video file; record =/=1 ---> don't record
        vidquality=100;  % quality of video, value between 0 to 100  (worst to best)
        pausetime=1; %number of seconds to pause at first frame in video
                % note:  ->for total number of frames (or time steps calculated), that is; totalframes=(time+pausetime)*fps
                %        -> Framerate of live output produced by matlab will not effect the
                %           framerate of the final produced video recording (which is purely defined by
                %           variable "fps"); so do not worry about the live matlab render
                %           being "choppy" or slow as the importance of this stage is to
                %           recorded each frame individually. 

                %           After the completion of the live render; the frames are put 
                %           together in a video file for which frames will be played at 
                %           rate as defined by variable "fps".

                %           Hence, you may need to run matlab for a long period to complete
                %           a record (i.e - a full or half of a day) depending on detail of
                %           output, some key parameters of influence are;
                %             --> "fps" ; increasing fps will increase total number of
                %                 frames(time steps) to calculate and record.
                %             --> "timelength"; same as above
                %             --> "numberofparticles"; increasing number of particles will
                %                 increase the time it takes to render each frame.

                %         -> if you only want to view the raw matlab output (and not
                %            interested in recording  ) set record to 0, as this will
                %            speed up the frame rate of the raw output.
        pauseb4record=0; % 1=following calculations, Matlab will wait for keyinput to begin video recording
                         %other=matlab automatically begins recording, following completion of calculations
        vidtitle='Myfirstmatlabvideo_gravity.avi'; % title of recorded video; 
         %note that if a file already exists with same name as vidtitle, this file will
         %be automatically over-written when the code is compiled
         diagramtitle='Rendered in Matlab with code created by Jamie M. Johns';


        % Further note: If the matlab render is interrupted (i.e-you stop the code
        % from completey rendering); you can save frames, that have already been
        % rendered, as a video file (with filename defined by vidtitle variable) by
        % running the following code in matlab window close(writerObj);


%###################################################################################################################################################


% CALCULATIONS BELOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        partnum=round(partnum);

        % setup figure$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            figure(1)
            if sum(outputpossize)~=0
            set(figure(1),'Position',outputpossize);
            end
            whitebg(figure(1), [0 0 0])
            if record==1
            writerObj = VideoWriter(vidtitle); %title
                 writerObj.FrameRate = fps; %frame rate
              writerObj.Quality = vidquality;   %100 equal 100 percent quality 
              open(writerObj);
            set(gca,'nextplot','replacechildren');
            set(gcf,'Renderer','zbuffer');
            end
        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        %window dialogue!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        commandwindow
        fprintf('Code executed at %.3f seconds\n\n',toc)
          fprintf('Details of parameters!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
          str1='Maxmass=%.4e | partnum=%.0f | G=%.4e\n';
           str2=' domainwidth=%.4e | elasticcoll=%.0f | ECP=%.4e\n'; 
           str3=' mindist=%.4e | fps=%.2f | timelength=%.2f\n';
           strg=[str1 str2 str3];
          fprintf(strg,Maxmass,partnum,G,domainwidth,elasticcoll,...
              ECP,mindist,fps,timelength)
          if record==1
           fprintf('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n ')
           fprintf('Output will be recorded with settings:\n ')
           fprintf('Output path-->%s \n',pwd)
           fprintf('filename:%s \n',vidtitle)
           pos = get(figure(1), 'Position'); 
           fprintf('resolution:%.0f by %.0f \n',pos(3),pos(4))
           fprintf('timelength:%0.f \n',timelength)
           fprintf('Frames per second:%0.f \n',fps)
           fprintf('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n ')
          end
        fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n')
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        commandwindow

        rotperframe=rotrate3d*360/(fps*60); %rate of rotating frame in 3d simulation
        dt=1/fps; %delta time

        x=rand(1,partnum,1)'.*domainwidth; %set random initial x position for each particle
        y=rand(1,partnum,1)'.*domainwidth; %set random initial y position for each particle
        z=rand(1,partnum,1)'.*domainwidth; %set random initial z position for each particle

        vel=zeros(partnum,3); %initial velocity (x,y,z) of each particle
        
        mass=Minmass+abs(Maxmass-Minmass).*rand(1,partnum); % generate random mass for each particle between min and max specified mass
        Mm=max(mass); %record maximum mass (used for graphics (plotting) to scale "size" of particles in visual output)




        %anonymous functions /////////////////////////////////////////////////////
        % overall used equation of force
        dist=@(A,B,C,D,E,F) (sqrt((A-B).^2+(C-D).^2+(E-F).^2));
        dist2=@(A,B,C,D,E,F,M) (dist(A,B,C,D,E,F)>=mindist).*dist(A,B,C,D,E,F)+(dist(A,B,C,D,E,F)<mindist).*mindist.*M;
        dist4=@(A,B,M) (abs(B-A)>=mindist).*(B-A)+(abs(B-A)<mindist).*mindist.*M;
        dist3=@(A,B,M) (abs(B-A)>=mindist).*abs(B-A)+(abs(B-A)<=mindist).*mindist.*M;
        vec=@(A,B,C,D,E,F,M) (5./(dist2(A,B,C,D,E,F,M).^2)).*((dist4(A,B,M))./dist3(A,B,M));
        F=@(A,B,C,D,E,F,n,l,M) ((G.*n.*l)./(1)).*vec(A,B,C,D,E,F,M); 

% DESCRIPTION OF  using equations of motion (above) and how calculations are performed
%                [description is brief due to me being very bust at the time,
        
         %F=G*m1*m2.*vec(x1,x2,y1,y2,z1,z2,C)  <---for x component of force particle 2 act on particle 1  [C is set to one by default]
          % vec(x1,x2,y1,y2,z1,z2)=dist4./(dist2*dist2*dist3)
          %             dist4=x2-x1  dist3=abs(x1-x2);
          %             dis2=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)

          %->For N particle/body system; j=1:N k=1:N       
          %--> for x-component of motion for the jth particle
             %Fxjk=G*mj*mk.*vec(xj,xk,yj,yk,zj,zk,C)
             %Fxj=sum[G*mj*mk.*vec(xj,xk,yj,yk,zj,zk,C); for k=1:N & k~=j]  [sum of all forces acting on jth particle]
             % also Fxj=mass(j)*accel(j) [newtons second law]
             %  and     --> accel(j)=(Vxj_final(j)-Vxj_initial(j))/(dt)     [Vxj = x-component velocity of jth particle]
             % therefore --->  Vxj_final(j)=Vxj_initial(j)+[{dt/mass(j)}*Fxj]
          %--> above is repeated for y and z component of jth particle
          %    and then repeated for x,y,z componentents of (j+1)th particle
          %-> above is completed for j=1:N particles for each time step
        dt=1/30;
        ind=0;
        tic



            [x1]=meshgrid(x(:,1)');
            [y1]=meshgrid(y(:,1)');
            [z1]=meshgrid(z(:,1)');
            [ma1,ma2]=meshgrid(mass,mass);
            onez=ones(size(ma1));
            Fx=zeros(partnum,partnum);
            Fy=zeros(partnum,partnum);
            Fz=zeros(partnum,partnum);

                 n=1:partnum;
                    m=1:partnum;



             fprintf('Time-step calculations started at %.3f seconds : \n',toc)
        for k=2:(pausetime+timelength).*fps;

            if k >= ((pausetime+timelength).*fps).*0.1 & ind==0;
            fprintf('%.0f%% of calculations completed at %.3f seconds \n',100.*k/((pausetime+timelength).*fps),toc)
            ind=1;
            elseif k >= ((pausetime+timelength).*fps).*0.2 & ind==1;
            fprintf('%.0f%% of calculations completed at %.3f seconds \n',100.*k/((pausetime+timelength).*fps),toc)
            ind=2;
                elseif k >= ((pausetime+timelength).*fps).*0.3 & ind==2;
            fprintf('%.0f%% of calculations completed at %.3f seconds \n',100.*k/((pausetime+timelength).*fps),toc)
            ind=3;
                   elseif k >= ((pausetime+timelength).*fps).*0.4 & ind==3;
            fprintf('%.0f%% of calculations completed at %.3f seconds \n',100.*k/((pausetime+timelength).*fps),toc)
            ind=ind+1; 
                   elseif k >= ((pausetime+timelength).*fps).*0.5 & ind==4;
            fprintf('%.0f%% of calculations completed at %.3f seconds \n',100.*k/((pausetime+timelength).*fps),toc)
            ind=ind+1; 
                       elseif k >= ((pausetime+timelength).*fps).*0.6 & ind==5;
            fprintf('%.0f%% of calculations completed at %.3f seconds \n',100.*k/((pausetime+timelength).*fps),toc)
            ind=ind+1; 
                       elseif k >= ((pausetime+timelength).*fps).*0.7 & ind==6;
            fprintf('%.0f%% of calculations completed at %.3f seconds \n',100.*k/((pausetime+timelength).*fps),toc)
            ind=ind+1; 
                       elseif k >= ((pausetime+timelength).*fps).*0.8 & ind==7;
            fprintf('%.0f%% of calculations completed at %.3f seconds \n',100.*k/((pausetime+timelength).*fps),toc)
            ind=ind+1; 
                       elseif k >= ((pausetime+timelength).*fps).*0.9 & ind==8;
            fprintf('%.0f%% of calculations completed at %.3f seconds \n',100.*k/((pausetime+timelength).*fps),toc)
            ind=ind+1; 
                       elseif k >= ((pausetime+timelength).*fps).*1 & ind==9;
            fprintf('%.0f%% of calculations completed at %.3f seconds \n',100.*k/((pausetime+timelength).*fps),toc)
            ind=ind+1; 

            else
            end



                x1=repmat(x(:,k-1)',partnum,1);
                y1=repmat(y(:,k-1)',partnum,1);
                z1=repmat(z(:,k-1)',partnum,1);
                x2=x1';
                y2=y1';
                z2=z1';

                    sclr=(k>=fps.*pausetime).*(pausetime~=0)+(pausetime==0);



                 Fx=Fx+(x1~=x2).*F(x1,x2,y2,y1,z2,z1,ma1,ma2,onez).*sclr;
                 Fy=Fy+(y1~=y2).*F(y1,y2,x2,x1,z2,z1,ma1,ma2,onez).*sclr;
                 Fz=Fz+(z1~=z2).*F(z1,z2,y2,y1,x2,x1,ma1,ma2,onez).*sclr;


                    Fy(isnan(Fy))=0;
                    Fx(isnan(Fx))=0;
                    Fz(isnan(Fx))=0; 
                Fz=(threedim==1).*Fz;



                 sfx(n,1)=sum(Fx(n,:));
                  sfy(n,1)=sum(Fy(n,:));   
                  sfz(n,1)=sum(Fz(n,:));  



         vel(n,3)=((sfz(n,1).*dt)./(mass(1,n)'))+vel(n,3);         
        vel(n,2)=((sfy(n,1).*dt)./(mass(1,n)'))+vel(n,2);
        vel(n,1)=((sfx(n,1).*dt)./(mass(1,n)'))+vel(n,1);


        x(n,k)=x(n,k-1)+vel(n,1).*dt;
        y(n,k)=y(n,k-1)+vel(n,2).*dt;
        z(n,k)=z(n,k-1)+vel(n,3).*dt;

        X=repmat(x(n,k),1,size(x(n,k),1));
        Y=repmat(y(n,k),1,size(x(n,k),1));
        Z=repmat(z(n,k),1,size(z(n,k),1));


        %if particle position exceeds boundary
             Fx=(Z<domainwidth).*(Z>0).*(Y<domainwidth).*(Y>0).*(X<domainwidth).*(X>0).*Fx;   
             Fy=(Z<domainwidth).*(Z>0).*(Y<domainwidth).*(Y>0).*(X<domainwidth).*(X>0).*Fy; 
             Fz=(Z<domainwidth).*(Z>0).*(Y<domainwidth).*(Y>0).*(X<domainwidth).*(X>0).*Fz; 
             if elasticcoll==1
                       vel(n,3)=(z(n,k)<domainwidth).*(z(n,k)>0).*vel(n,3)-((z(n,k)>=domainwidth)+(z(n,k)<=0)).*vel(n,3).*ECP;  
              vel(n,2)=(y(n,k)<domainwidth).*(y(n,k)>0).*vel(n,2)-((y(n,k)>=domainwidth)+(y(n,k)<=0)).*vel(n,2).*ECP;  
              vel(n,1)=(x(n,k)<domainwidth).*(x(n,k)>0).*vel(n,1)-((x(n,k)>=domainwidth)+(x(n,k)<=0)).*vel(n,1).*ECP;
             else
            vel(n,1)=(z(n,k)<domainwidth).*(z(n,k)>0).*(y(n,k)<domainwidth).*(y(n,k)>0).*(x(n,k)<domainwidth).*(x(n,k)>0).*vel(n,1);   
            vel(n,2)=(z(n,k)<domainwidth).*(z(n,k)>0).*(x(n,k)<domainwidth).*(x(n,k)>0).*(y(n,k)<domainwidth).*(y(n,k)>0).*vel(n,2);    
               vel(n,3)=(z(n,k)<domainwidth).*(z(n,k)>0).*(x(n,k)<domainwidth).*(x(n,k)>0).*(y(n,k)<domainwidth).*(y(n,k)>0).*vel(n,3);    
             end
        x(n,k)=(x(n,k)>=domainwidth).*domainwidth.*0.999+(x(n,k)<0).*0.0001+(x(n,k)<=domainwidth).*(x(n,k)>0).*x(n,k);
        y(n,k)=(y(n,k)>=domainwidth).*domainwidth.*0.999+(y(n,k)<0).*0.0001+(y(n,k)<=domainwidth).*(y(n,k)>0).*y(n,k);   
        z(n,k)=(z(n,k)>=domainwidth).*domainwidth.*0.999+(z(n,k)<0).*0.0001+(z(n,k)<=domainwidth).*(z(n,k)>0).*z(n,k); 
        end 
        fprintf('calculations completed at %.3f seconds \n',toc)
        if record==1
           if pauseb4record==1 
            WAITING=input('Type in any key and press enter to start render and recording:');  
           end
        end    
        ind=0
        xlim([0 domainwidth])
        ylim([0 domainwidth])
        hold on

        fprintf('Inital plotting started at %.3f seconds \n',toc)

        if showdetail==1 && record==1;
            detailtime=(detailtime<=1)+(detailtime>1)*(detailtime);
          text(domainwidth.*(3/10),domainwidth.*(9/10),'Details of parameters:','fontsize',25,'color',[0 1 0])
          str1='Maxmass=%.4e | partnum=%.0f | G=%.4e\n';
           str2=' domainwidth=%.4e | elasticcoll=%.0f | ECP=%.4e\n'; 
           str3=' mindist=%.4e | fps=%.2f | timelength=%.2f\n';
           strg=[str1 str2 str3];
          text(domainwidth.*(1/10),domainwidth.*(7.5/10),sprintf(strg,Maxmass,partnum,G,domainwidth,elasticcoll,...
              ECP,mindist,fps,timelength),'fontsize',15)

         for k=1:detailtime*fps 
             drawnow
                   frame = getframe(figure(1));
                 writeVideo(writerObj,frame);    

         end   
        end


        clf(1)
        xlim([0 domainwidth])
        ylim([0 domainwidth])
        hold on


        if threedim==1
         for k=randperm(partnum)
           pl(1,k)=plot3(x(k,2),y(k,2),z(k,2),'r.','Markers',25+45.*(mass(1,k)./Mm),'color',[mass(1,k)/Mm 0.5 rand(1)]);  
        end       
        else
        for k=randperm(partnum)
           pl(1,k)=plot(x(k,2),y(k,2),'r.','Markers',25+45.*(mass(1,k)./Mm),'color',[mass(1,k)/Mm 0.5 rand(1)]);  
        end 
        end
        drawnow
        view(2)



        fprintf('Inital plotting complete at and simulation plot started at %.3f seconds \n',toc)






        fprintf('\n Simulation plot started at %.3f seconds \n',toc)

        RT=0;
        for t=1:(pausetime+timelength).*fps;
                if t > ((pausetime+timelength).*fps).*0.1 & ind==0;
            fprintf('%.0f%% of plotting completed at %.3f seconds \n',100.*t/((pausetime+timelength).*fps),toc)
            ind=1;
            elseif t > ((pausetime+timelength).*fps).*0.2 & ind==1;
            fprintf('%.0f%% of plotting completed at %.3f seconds \n',100.*t/((pausetime+timelength).*fps),toc)
            ind=2;
                elseif t > ((pausetime+timelength).*fps).*0.3 & ind==2;
            fprintf('%.0f%% of plotting completed at %.3f seconds \n',100.*t/((pausetime+timelength).*fps),toc)
            ind=3;
                   elseif t > ((pausetime+timelength).*fps).*0.5 & ind==3;
            fprintf('%.0f%% of plotting completed at %.3f seconds \n',100.*t/((pausetime+timelength).*fps),toc)
            ind=4; 
                       elseif t > ((pausetime+timelength).*fps).*0.75 & ind==4;
            fprintf('%.0f%% of plotting completed at %.3f seconds \n',100.*t/((pausetime+timelength).*fps),toc)
            ind=5; 
                    elseif t > ((pausetime+timelength).*fps).*0.9 & ind==5;
            fprintf('%.0f%% of plotting completed at %.3f seconds \n',100.*t/((pausetime+timelength).*fps),toc)
            ind=6;  
                     elseif t == ((pausetime+timelength).*fps)& ind==6;
            fprintf('%.0f%% of plotting completed at %.3f seconds \n',100.*t/((pausetime+timelength).*fps),toc)
            ind=7;  
            else
                end
             if threedim==1
         for k=1:partnum

           set(pl(1,k),'Xdata',x(k,t),'Ydata',y(k,t),'zdata',z(k,t))
         end
         view(RT,22)
         axis([0 domainwidth 0 domainwidth 0 domainwidth])
             else
                  for k=1:partnum
                    set(pl(1,k),'Xdata',x(k,t),'Ydata',y(k,t)) 
                  end
                     view(2)
                      axis([0 domainwidth 0 domainwidth])
             end 

        title(diagramtitle,'fontsize',15)
          drawnow  
              if record==1
                 frame = getframe(figure(1));
                 writeVideo(writerObj,frame);
              end
                 RT=RT+rotperframe;
        end 

        fprintf('Output completed at %.3f seconds \n',toc)
        if record==1
        close(writerObj);
        end



