function [Particle] = MPM_solver(SolidModel, Cell, Node, Particle, Time, Physics)

%% Calculate the value of shape function and gradient of the shape function
x_test = Particle.x;
[Particle,Cell] = Compute_Interpolator_MPM(Particle,Cell,Node,x_test); 

 %% Mapping from particle to nodes
[Node]=Interpolate_Particle_To_Grid(Node,Particle);

%% Update momentum and acceleration
[Node]=Integration(Node,Time,Physics,Particle);

%% Update solid particle velocity and position
[Particle] = Update_Particle(Particle,Node,Time);

%% Update constitutive model
 switch SolidModel.name
            case 'Neo_Hookean_Elastic'
                [Particle]=Neo_Hookean_elastic(SolidModel,Particle);            
 end        
