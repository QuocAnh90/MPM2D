close all; clear all;
% Unit: Newton - seconds - metre

addpath('SubFunctions');
addpath('SubFunctions/Cores');
addpath('SubFunctions/Constitutive_Models');
addpath('SubFunctions/CPDI_solvers');
addpath('SubFunctions/MPM_solvers');
addpath('SubFunctions/GIMP_solvers');

%% Please select the versions of MPM!!!!!!!!!!!!!!!!!
% Original MPM: 'MPM' or 'GIMP' or 'CPDI'
Version = 'GIMP';

%% Constutitive model input
% Neo_Hookean_Elastic % Mohr_Coulomb & Mohr_Coulomb_rotation
SolidModel.name = 'Mohr_Coulomb_rotation';
Physics.gravity             = [0 -10]           ;                   % gravity acceleration
SolidModel.density          = 1500.0            ;                   % solid density

% Material porperties for Neo_Hookean_Elastic
% SolidModel.Young_modul      = 1000000           ;                   % Young modulus of solid
% SolidModel.nu               = 0.3               ;                   % Poison ratio

%% Material porperties for Mohr_Coulomb
SolidModel.Young_modul      = 10000000           ;                   % Young modulus of solid
SolidModel.nu               = 0.3               ;                   % Poison ratio
SolidModel.Friction         = 30                ;                   % Friction angle
SolidModel.Dilation         = 0.0               ;                   % Dilation angle
SolidModel.cohesion         = 0.3               ;                   % cohesion

%% Structured Grid input
% Node.CountX                = 67;                % number of nodes in X direction
% Node.CountY                = 27;                % number of nodes in Y direction
% Cell.size(1)               = 0.1;               % size of element in X direction
% Cell.size(2)               = 0.1;               % size of element in Y direction
% % Boundary coordination (symmetry)
% x_min = 0.1; x_max = 6.5; y_min = 0.2; y_max = 2.5;

Node.CountX                = 134;                % number of nodes in X direction
Node.CountY                = 54;                % number of nodes in Y direction
Cell.size(1)               = 0.05;               % size of element in X direction
Cell.size(2)               = 0.05;               % size of element in Y direction
% Boundary coordination (symmetry)
x_min = 0.05; x_max = 6.5; y_min = 0.1; y_max = 2.5;

%% Time input
Time.wavespeed          = sqrt(SolidModel.Young_modul/SolidModel.density);
Time.finalTime          = 5.00;
Time.timestep           = 0.00005;
RealTime                = 0;

%% Particle input
Particle.PPCX      = 2;                                 % Particle Per Cell
Particle.PPCY      = 2;                                 % Particle Per Cell
Particle.PPC       = Particle.PPCX * Particle.PPCY;     % Particle Per Cell
Particle.Count     = 800*Particle.PPC;                   % Total number of particle

%% Grid generation
[Node,Cell] = Grid_Generation(Node,Cell,x_min,x_max,y_min,y_max);

%% Boundary condition
% for n=1:Node.Count
%     if Node.x(n,1)<=x_min 
%         Node.Count_BoundaryY = Node.Count_BoundaryY+1;
%         Node.BoundaryY = [Node.BoundaryY n];
%     end
% end

for n=1:Node.Count
    if Node.x(n,2)<=y_min 
        Node.Count_BoundaryX = Node.Count_BoundaryX+1;
        Node.BoundaryX = [Node.BoundaryX n];
    end
end

%% Find interacting nodes of each cell
[Cell,Particle] = ConnectNode(Cell, Node,Particle,Version);

%% Particle generation
% Create variables
[Particle] = Particle_Generation(Particle,Cell,Node,SolidModel);

% Generate particles
sp=1;
while sp<Particle.Count+0.0001
    for i=1:80
        for j=1:40
            Particle.x(sp,1:2)= [1*Cell.size(1)+0.5*Particle.size(1)+(j-1)*Particle.size(1) 2*Cell.size(2)+0.5*Particle.size(2)+(i-1)*Particle.size(2)];
            sp=sp+1;
        end
    end
end
Particle.x_ini          = Particle.x;                           % initial position
                 
%% Plot initial condition
initial_figure = Plot_Initial(Particle.x,Node.x,Cell.size);

%% start the algorithm
% video
step = 100;     % number of frame to save
r=step/5;      % number of frame per second video ~200s

writerObj2           = VideoWriter('column.avi');
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
        
        case 'GIMP'
        [Particle] = GIMP_solver(SolidModel, Cell, Node, Particle, Time, Physics);

         case 'CPDI'
        [Particle] = CPDI_solver(SolidModel, Cell, Node, Particle, Time, Physics);
      end
     
     % Update timestep 
     RealTime = RealTime+Time.timestep;
 end

 %% Generate video of the result
    Figure=Plot_Final(Particle,Node,Cell);
    frame2 = getframe(Figure);
    writeVideo(writerObj2,frame2);
    
    end
    
close(writerObj2);
  