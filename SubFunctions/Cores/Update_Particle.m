function [Particle] = Update_Particle(Particle,Node,Time)

%% Update particle velocity and position
AccelerationX   = Node.acceleration(:,1);
AccelerationY   = Node.acceleration(:,2);
VelocityX       = Node.velocity(:,1);
VelocityY       = Node.velocity(:,2);
dt              = Time.timestep;             % Time step

% Update particle's velocity
Particle.velocity(:,1)      = Particle.velocity(:,1) + sum(Particle.S.*AccelerationX(Particle.CONNECT),2)*dt;
Particle.velocity(:,2)      = Particle.velocity(:,2) + sum(Particle.S.*AccelerationY(Particle.CONNECT),2)*dt;

% Update particle's position
Particle.x(:,1)             = Particle.x(:,1) + sum(Particle.S.*VelocityX(Particle.CONNECT),2)*dt;
Particle.x(:,2)             = Particle.x(:,2) + sum(Particle.S.*VelocityY(Particle.CONNECT),2)*dt;

% Update particle's velocity gradient
Particle.Gradvelocity(:,1)  = sum(Particle.dSx.*VelocityX(Particle.CONNECT),2);
Particle.Gradvelocity(:,2)  = sum(Particle.dSx.*VelocityY(Particle.CONNECT),2);
Particle.Gradvelocity(:,3)  = sum(Particle.dSy.*VelocityX(Particle.CONNECT),2);
Particle.Gradvelocity(:,4)  = sum(Particle.dSy.*VelocityY(Particle.CONNECT),2);

% Update particle's deformation gradient
% F11 = dv11*F11 + dv12*F21
% F12 = dv11*F12 + dv12*F22
% F21 = dv21*F11 + dv22*F21
% F22 = dv21*F12 + dv22*F22
Particle.defgrad_Old = Particle.defgrad;
Particle.defgrad (:,1)  = (1+Particle.Gradvelocity(:,1)*dt).*Particle.defgrad_Old(:,1) + (Particle.Gradvelocity(:,2)*dt).*Particle.defgrad_Old(:,3);
Particle.defgrad (:,2)  = (1+Particle.Gradvelocity(:,1)*dt).*Particle.defgrad_Old(:,2) + (Particle.Gradvelocity(:,2)*dt).*Particle.defgrad_Old(:,4);
Particle.defgrad (:,3)  = (Particle.Gradvelocity(:,3)*dt).*Particle.defgrad_Old(:,1) + (1+Particle.Gradvelocity(:,4)*dt).*Particle.defgrad_Old(:,3);
Particle.defgrad (:,4)  = (Particle.Gradvelocity(:,3)*dt).*Particle.defgrad_Old(:,2) + (1+Particle.Gradvelocity(:,4)*dt).*Particle.defgrad_Old(:,4);

% Update particle's volume
Particle.Jacobian = Particle.defgrad (:,1).*Particle.defgrad (:,4) - Particle.defgrad (:,2).*Particle.defgrad (:,3);
Particle.volume = Particle.Jacobian .* Particle.volume_ini;