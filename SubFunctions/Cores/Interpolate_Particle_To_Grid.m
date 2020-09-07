function [Node]=Interpolate_Particle_To_Grid(Node,Particle)

%% Input
S = Particle.S;                 % Shape fucntion (size Np x Nn)
dSx = Particle.dSx;             % Shape function gradient in x direction 
dSy = Particle.dSy;             % Shape function gradient in y direction
mp = Particle.mass;             % Mass
vp = Particle.velocity;         % Velocity
Nn = Particle.Node;             % Number of interacting nodes
Np = Particle.Count;            % Total number of particles

%% Initialization of nodal variables
Node.mass(:)        = 0;         % Nodal Mass
Node.momentum(:)    = 0;         % Nodal Momentum
Node.inForce(:)     = 0;         % Nodal Internal force
% Node.exForce(:)     = 0;         % Nodal External force

%% Reshape particle variables 
mass        = reshape(S.*repmat(mp,1,Nn)            ,Np*Nn,1); 
momentumX   = reshape(S.*repmat(mp.*vp(:,1),1,Nn)   ,Np*Nn,1); 
momentumY   = reshape(S.*repmat(mp.*vp(:,2),1,Nn)   ,Np*Nn,1); 

%% Interpolation from particle to grid task
Node.nmass              = accumarray(Particle.CONNECT(:) ,mass   ,[Node.Count 1]);      
Node.nmomentum(:,1)     = accumarray(Particle.CONNECT(:) ,momentumX   ,[Node.Count 1]);      
Node.nmomentum(:,2)     = accumarray(Particle.CONNECT(:) ,momentumY   ,[Node.Count 1]);      

 