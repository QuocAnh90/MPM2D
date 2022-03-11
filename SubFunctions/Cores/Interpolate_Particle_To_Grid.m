function [Node]=Interpolate_Particle_To_Grid(Node,Particle)

%% Input
S = Particle.S;                 % Shape fucntion (size Np x Nn)
% dSx = Particle.dSx;             % Shape function gradient in x direction 
% dSy = Particle.dSy;             % Shape function gradient in y direction
mp = Particle.mass;             % Mass
vp = Particle.velocity;         % Velocity
Nnp = Particle.Node;            % Number of interacting nodes
Np = Particle.Count;            % Total number of particles
Nn = Node.Count;                % Total number of nodes
DOF = Node.DOF;                 % Number of dimension

%% Initialization of nodal variables
Node.mass(:)        = 0;         % Nodal Mass
Node.momentum(:)    = 0;         % Nodal Momentum
Node.inForce(:)     = 0;         % Nodal Internal force
% Node.exForce(:)     = 0;         % Nodal External force

%% Reshape particle variables 
mass        = reshape(S.*repmat(mp,1,Nnp)            ,Np*Nnp,1); 
% momentum    = reshape([S.*repmat(mp.*vp(:,1),1,Nnp) ; S.*repmat(mp.*vp(:,2),1,Nnp)], Np*Nnp*DOF,1);
momentumX   = reshape(S.*repmat(mp.*vp(:,1),1,Nnp)   ,Np*Nnp,1); 
momentumY   = reshape(S.*repmat(mp.*vp(:,2),1,Nnp)   ,Np*Nnp,1); 

inForce     = squeeze(sum(...
                        Particle.B .* repmat(...
                            reshape(...
                                Particle.stress,size(Particle.stress,1),1,Np),1,Nnp*DOF),1))...
              .* repmat(Particle.volume,1,Nnp*2);
        
%% Interpolation from particle to grid task
Node.mass               = accumarray(Particle.CONNECT(:) ,mass                      ,[Nn 1]);    

% CONNECT  = [DOF*Particle.CONNECT-1;DOF*Particle.CONNECT];
% Node.momentum           = accumarray(CONNECT             ,momentum                  ,[Nn*DOF 1]);   
% Node.momentum should be 2Nn x 1 instead of Nn Nn

Node.momentum(:,1)      = accumarray(Particle.CONNECT(:) ,momentumX                 ,[Nn 1]);      
Node.momentum(:,2)      = accumarray(Particle.CONNECT(:) ,momentumY                 ,[Nn 1]);  

for i =1:Nnp
Node.inForce(:,1)       = Node.inForce(:,1) - accumarray(Particle.CONNECT(:,i)   ,inForce(2*i-1,:)    ,[Nn 1]);
Node.inForce(:,2)       = Node.inForce(:,2) - accumarray(Particle.CONNECT(:,i)   ,inForce(2*i,:)      ,[Nn 1]);
end

test = 1;

%% Document the code
% For calculation of internal force
% Particle.stress (3 x Np) 
% Reshape paticle.stress (3 x 1 x Np)
% Rempat particle stress (3 x 8 x Np)
% sum ( B(3 x 8 x Np) .* stress, 1)  = 1 x 8 x Np and squeeze
