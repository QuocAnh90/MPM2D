function [Particle] = GIMP_solver(SolidModel, Cell, Node, Particle, Time, Physics)

%% Calculate the value of shape function and gradient of the shape function
[Particle,Cell] = Compute_Interpolator_GIMP(Particle,Cell,Node); 

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
            case 'Mohr_Coulomb'
                [Particle]=Mohr_Coulomb(SolidModel,Particle,Time);
            case 'Mohr_Coulomb_rotation'
                [Particle]=Mohr_Coulomb_rotation(SolidModel,Particle,Time);
 end        
