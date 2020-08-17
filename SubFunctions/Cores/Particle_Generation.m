function [Particle] = Particle_Generation(Particle,Cell,SolidModel,gravity)

Count = Particle.Count;

%% Particle variables
Particle.size(1)        = Cell.size(1)/sqrt(Particle.PPC);      % size of particle in X direction
Particle.size(2)        = Cell.size(2)/sqrt(Particle.PPC);      % size of particle in Y direction
Particle.x              = zeros(Count,2);                       % Position
Particle.d              = zeros(Count,2);                       % Displacement
Particle.x_ini          = Particle.x;                           % initial position
Particle.density        = SolidModel.density * ones(Count,1);   % Density
Particle.body           = [zeros(Count,1) -gravity*ones(Count,1)];    % body force
Particle.stress         = zeros(Count,3);                       % stress
Particle.velocity       = zeros(Count,2);                       % velocity
Particle.strain         = zeros(Count,3);                       % strain
Particle.traction       = zeros(Count,2);                       % traction
Particle.defgrad        = cell(Count,2);                        % Deformation gradient
Particle.volume         = zeros(Count,1);                       % Volume
Particle.r1             = zeros(Count,2);
Particle.r2             = zeros(Count,2);

%% Initial condition
% Initialize the gradient deformation and vector r1 and r2 (for CPDI)
for p = 1:Count
    Particle.r1(p,:) = [Particle.size(1)/2 0];
    Particle.r2(p,:) = [0 Particle.size(2)/2];
    Particle.defgrad{p} = [1 0; 0 1];
end
Particle.r1_ini = Particle.r1;
Particle.r2_ini = Particle.r1;

% Compute the particle volume based on r1 and r2
for p=1:Count
Particle.volume(p)      = 4*abs(Particle.r1(p,1)*Particle.r2(p,2)-Particle.r1(p,2)*Particle.r2(p,1)); 
end
Particle.volume_ini     = Particle.volume;
Particle.mass           = Particle.density .* Particle.volume;   % mass
