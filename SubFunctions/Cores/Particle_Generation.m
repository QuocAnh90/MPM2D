function [Particle] = Particle_Generation(Particle,Cell,Node,SolidModel)

Count = Particle.Count;
DOF = Node.DOF;                 % Number of dimension

%% Particle variables
Particle.size(1)        = Cell.size(1)/Particle.PPCX;      % size of particle in X direction
Particle.size(2)        = Cell.size(2)/Particle.PPCY;      % size of particle in Y direction

% Scalar
Particle.volume         = zeros(Count,1);                       % Volume

% Vector
Particle.x              = zeros(Count,2);                       % Position
Particle.displacement   = zeros(Count,2);                       % Displacement
Particle.density        = SolidModel.density * ones(Count,1);   % Density
Particle.velocity       = zeros(Count,2);                       % velocity
Particle.Gradvelocity   = zeros(Count,4);                       % Gradient velocity
% Particle.traction       = zeros(Count,2);                     % traction
Particle.defgrad        = zeros(Count,4);                       % Deformation gradient  (11, 12, 21, 22)
Particle.r1             = zeros(Count,2);
Particle.r2             = zeros(Count,2);

% Tensor
Particle.strainInc      = zeros(Count,3);                     % strain increment (xx, yy, xy)
Particle.stress         = zeros(3,Count);                       % stress (xx, yy, xy)
Particle.B              = zeros(3,DOF * Particle.Node,Count);   % B matrix

%% Initial condition
% Initialize the gradient deformation and vector r1 and r2 (for CPDI)
for p = 1:Count
    Particle.r1(p,:) = [Particle.size(1)/2 0];
    Particle.r2(p,:) = [0 Particle.size(2)/2];
end
Particle.defgrad = [ones(Count,1) zeros(Count,1) zeros(Count,1) ones(Count,1)];
Particle.defgrad_Old = Particle.defgrad;
Particle.r1_ini = Particle.r1;
Particle.r2_ini = Particle.r1;

% Compute the particle volume based on r1 and r2
for p=1:Count
Particle.volume(p)      = 4*abs(Particle.r1(p,1)*Particle.r2(p,2)-Particle.r1(p,2)*Particle.r2(p,1)); 
end
Particle.volume_ini     = Particle.volume;
Particle.mass           = Particle.density .* Particle.volume;   % mass
