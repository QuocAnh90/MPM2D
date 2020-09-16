close all; clear all;
% Unit: Newton - seconds - metre

%% Please select the versions of MPM!!!!!!!!!!!!!!!!!
% Original MPM: 'MPM'; % CPDI: 'CPDI'
Version = 'MPM';

%% Constutitive model input
% Linear_Elastic % Neo_Hookean_Elastic
SolidModel.name = 'Neo_Hookean_Elastic';

% Material porperties
SolidModel.Young_modul      = 1000000           ;                   % Young modulus of solid
SolidModel.density          = 1050.0            ;                   % solid density
SolidModel.nu               = 0.3               ;                   % Poison ratio
Physics.gravity             = [0 -10]           ;                   % gravity acceleration

%% Structured Grid input
Node.CountX                = 6;                % number of nodes in X direction
Node.CountY                = 5;                % number of nodes in Y direction
Cell.size(1)               = 0.1;               % size of element in X direction
Cell.size(2)               = 0.1;               % size of element in Y direction
% Boundary coordination
x_min = 0; x_max = 0.5; y_min = 0; y_max = 0.4;

%% Time input
Time.wavespeed          = sqrt(SolidModel.Young_modul/SolidModel.density);
Time.finalTime          = 5;
Time.timestep           = 0.001;
RealTime                = 0;

%% Particle input
Particle.PPC       = 4;                 % Particle Per Cell
Particle.Count     = 2*Particle.PPC;   % Total number of particle

%% Grid generation
[Node,Cell] = Grid_Generation(Node,Cell,x_min,x_max,y_min,y_max);

%% Find interacting nodes of each cell
[Cell,Particle] = ConnectNode(Cell, Node,Particle);

%% Particle generation
% Create variables
[Particle] = Particle_Generation(Particle,Cell,SolidModel,Physics.gravity);

% Generate particles
sp=1;
while sp<Particle.Count+0.0001
    for i=1:4
        for j=1:2
            Particle.x(sp,1:2)= [2*Cell.size(1)+0.5*Particle.size(1)+(j-1)*Particle.size(1) 1*Cell.size(2)+0.5*Particle.size(2)+(i-1)*Particle.size(2)];
            sp=sp+1;
        end
    end
end
Particle.velocity = [1 2; 3 4; 5 6; 7 8;
                     1 2; 3 4; 5 6; 7 8];

Particle.stress   = [1 2 3 4; 5 6 7 8;
                     1 2 3 4; 5 6 7 8;
                     1 2 3 4; 5 6 7 8;
                     1 2 3 4; 5 6 7 8;];
                 
%% Plot initial condition
initial_figure = Plot_Initial(Particle.x,Node.x,Cell.size);

%% start the algorithm
% video
step = 100;     % number of frame to save
r=step/20;      % number of frame per second video ~200s

writerObj2           = VideoWriter('BendingBeam.avi');
writerObj2.FrameRate = r;    % number of frame per second
open(writerObj2);

    for tt = 1:step
    ft              = Time.finalTime/step*tt;
%     ft=Time.finalTime;

 while RealTime<ft+0.0000000001      
     RealTime
     
     switch Version
         case 'MPM'
        [Particle] = MPM_solver(SolidModel, Cell, Node, Particle, Time, Physics);

         case 'CPDI'
        [Particle] = CPDI_solver(SolidModel, Cell, Node, Particle, Time, Physics);
      end
     
     % Update timestep 
     RealTime = RealTime+dt;
 end

 %% Generate video of the result
    Figure=Plot_Final(Particle,Node,Cell);
    frame2 = getframe(Figure);
    writeVideo(writerObj2,frame2);
    
    end
    
close(writerObj2);
  