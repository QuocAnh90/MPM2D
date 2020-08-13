close all;clear all;
% Unit
% Newton - seconds - metre

%% Please select the versions of MPM!!!!!!!!!!!!!!!!!
% Original MPM: 'MPM'; % CPDI: 'CPDI'
Version = 'MPM';

%% Constutitive model
% Linear_Elastic % Neo_Hookean_Elastic
CM.name = 'Neo_Hookean_Elastic';

% Material porperties
CM.Young_modul      = 1000000           ;                   % Young modulus of solid
CM.psp              = 1050.0            ;                   % solid density
CM.nu               = 0.3               ;                   % Poison ratio
gravity                 = 10.0          ;                   % gravity acceleration

%% Structured Grid input
Node.NN(1)                 = 13;                % number of nodes in X direction
Node.NN(2)                 = 13;                % number of nodes in Y direction
Cell.le(1)                 = 0.5;               % size of element in X direction
Cell.le(2)                 = 0.5;               % size of element in Y direction
% Boundary coordination
x_min = 1; x_max = 6; y_min = 0; y_max = 6;

%% Grid generation
[Node,Cell] = Grid_Generation(Node,Cell,x_min,x_max,y_min,y_max);

%% Time
Time.wavespeed          = sqrt(CM.Young_modul/CM.psp);
Time.finalTime          = 5;
Time.timestep           = 0.001;
RealTime                = 0;

%% Particle input
Particle.PPC       = 9;                 % Particle Per Cell
Particle.Count     = 16*Particle.PPC;   % Total number of particle

%% Particle generation
[Particle] = Particle_Generation(Particle);


spCount                 = 16*PPC;
lp(1)                   = le(1)/sqrt(PPC);                                 % size of particle in X direction
lp(2)                   = lp(1);                              % size of particle in Y direction
x_sp                    = zeros(spCount,1);
d_sp                    = zeros(spCount,2);

sp=1;
while sp<spCount+0.0001
    for i=1:6
        for j=1:24
            x_sp(sp,1:2)= [2*le(1,1)+0.5*lp(1,1)+(j-1)*lp(1,1) 8*le(1,2)+0.5*lp(1,2)+(i-1)*lp(1,2)];
            sp=sp+1;
        end
    end
end
% x_sp: Vector, position of MPs
% spCount: total number of MPs

%% Plot initial condition
initial_figure = Plot_Initial(x_sp,LOC,le);

%% Particle variables
dparticle               = lp(1)*lp(2);                          % area of particle domain
x_spo                   = x_sp;                                 % initial position
p_sp                    = psp * ones(spCount,1);                % Density
b_sp                    = [zeros(spCount,1) -g*ones(spCount,1)];% body force

d_sp                    = zeros(spCount,2);                     % displacement
s_sp                    = zeros(spCount,3);                     % Stress tensor
ds_sp                   = zeros(spCount,3);                     % Stress increment
v_ssp                   = zeros(spCount,2);                     % velocty
e_sp                    = zeros(spCount,3);                     % Strain tensor
de_sp                   = zeros(spCount,3);                     % Strain increment
ptraction_sp            = zeros(spCount,2);                     % traction
F_sp                    = cell(spCount,1);                      % Gradient deformation
r1_sp                   = zeros(spCount,2);
r2_sp                   = zeros(spCount,2);

%% Initial condition
% Gradient deformation
for sp = 1:spCount
    r1_sp(sp,:) = [lp(1,1)/2 0];
    r2_sp(sp,:) = [0 lp(1,2)/2];
    F_sp{sp} = [1 0; 0 1];
end
r10_sp = r1_sp;
r20_sp = r2_sp;
V_sp                    = zeros(spCount,1);
for sp=1:spCount
V_sp(sp)                = 4*abs(r1_sp(sp,1)*r2_sp(sp,2)-r1_sp(sp,2)*r2_sp(sp,1)); 
end
V_spo                   = V_sp;
m_sp                    = psp * V_sp;                           % mass

% Traction
for sp=1:spCount
    if x_sp(sp,1)>1.01
        ptraction_sp(sp,1) = 1;
    end
end

%% start the algorithm
% video
step = 100;     % number of frame to save
r=step/20;      % number of frame per second video ~200s

writerObj2           = VideoWriter('BendingBeam.avi');
writerObj2.FrameRate = r;    % number of frame per second
open(writerObj2);

    for tt = 1:step
    Time.finalTime              = ftime/step*tt;
%     Time.finalTime=ftime;
 while RealTime<ft+0.0000000001      
     RealTime
     
     switch Version
         case 'MPM'
        [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp] = MPM_solver(CM,CModel_parameter,...
    nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
    nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt);

         case 'CPDI'
        [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,r1_sp,r2_sp] = CPDI_solver(CM,CModel_parameter,...
            nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
            nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt,r1_sp,r10_sp,r2_sp,r20_sp);
      end
     
     
 % Update time and step 
 RealTime = RealTime+dt;
 end

 %% Plot the result
StressProfile1=Plot_Final(x_sp,LOC,le,d_sp,e_sp,v_ssp,spCount,r1_sp,r2_sp);

    frame2 = getframe(StressProfile1);
    writeVideo(writerObj2,frame2);
    end
    
    close(writerObj2);
  