function [Particle] = Update_Particle(Particle,Node,Time)

%% Update particle velocity and position
AccelerationX   = Node.acceleration(:,1);
AccelerationY   = Node.acceleration(:,2);
VelocityX       = Node.velocity(:,1);
VelocityY       = Node.velocity(:,2);
dt              = Time.timestep;             % Time step

% UPdate particle's velocity
Particle.velocity(:,1)      = Particle.velocity(:,1) + sum(Particle.S.*AccelerationX(Particle.CONNECT),2)*dt;
Particle.velocity(:,2)      = Particle.velocity(:,2) + sum(Particle.S.*AccelerationY(Particle.CONNECT),2)*dt;

%Update particle's position
Particle.x(:,1)             = Particle.x(:,1) + sum(Particle.S.*VelocityX(Particle.CONNECT),2)*dt;
Particle.x(:,2)             = Particle.x(:,2) + sum(Particle.S.*VelocityY(Particle.CONNECT),2)*dt;
