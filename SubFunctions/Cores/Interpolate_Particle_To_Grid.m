function [Node]=Interpolate_Particle_To_Grid(Node,Particle)

%% Input
S = Particle.S;                 % Shape fucntion (size Np x Nn)
dSx = Particle.dSx;             % Shape function gradient in x direction 
dSy = Particle.dSy;             % Shape function gradient in y direction
mp = Particle.mass;             % Mass
vp = Particle.velocity;         % Velocity
Nnp = Particle.Node;            % Number of interacting nodes
Np = Particle.Count;            % Total number of particles
Nn = Node.Count;                % Total number of nodes

%% Initialization of nodal variables
Node.mass(:)        = 0;         % Nodal Mass
Node.momentum(:)    = 0;         % Nodal Momentum
Node.inForce(:)     = 0;         % Nodal Internal force
% Node.exForce(:)     = 0;         % Nodal External force

%% Reshape particle variables 
mass        = reshape(S.*repmat(mp,1,Nnp)            ,Np*Nnp,1); 
momentumX   = reshape(S.*repmat(mp.*vp(:,1),1,Nnp)   ,Np*Nnp,1); 
momentumY   = reshape(S.*repmat(mp.*vp(:,2),1,Nnp)   ,Np*Nnp,1); 
inForce     = squeeze(sum(Particle.B.*...
            repmat(permute(Particle.stress,[3 1 2]),Nnp*2,1),3).*...
            repmat(permute(Particle.volume,[3 1 2]),Nnp*2,1));
        
%% Interpolation from particle to grid task
Node.mass               = accumarray(Particle.CONNECT(:) ,mass                      ,[Nn 1]);      
Node.momentum(:,1)      = accumarray(Particle.CONNECT(:) ,momentumX                 ,[Nn 1]);      
Node.momentum(:,2)      = accumarray(Particle.CONNECT(:) ,momentumY                 ,[Nn 1]);  

for i =1:Nnp
Node.inForce(:,1)       = Node.inForce(:,1) - accumarray(Particle.CONNECT(:,i)   ,inForce(2*i-1,:)    ,[Nn 1]);
Node.inForce(:,2)       = Node.inForce(:,2) - accumarray(Particle.CONNECT(:,i)   ,inForce(2*i,:)      ,[Nn 1]);
end
