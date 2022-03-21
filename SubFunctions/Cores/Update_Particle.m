function [Particle] = Update_Particle(Particle,Node,Time)

%% Update particle velocity and position
AccelerationX   = Node.acceleration(:,1);
AccelerationY   = Node.acceleration(:,2);
VelocityX       = Node.velocity(:,1);
VelocityY       = Node.velocity(:,2);
dt              = Time.timestep;             % Time step

% Update particle's velocity
Particle.velocity(:,1)      = Particle.velocity(:,1) + sum(Particle.S.*AccelerationX(Particle.CONNECT),2)*dt;   % particle velocity x
Particle.velocity(:,2)      = Particle.velocity(:,2) + sum(Particle.S.*AccelerationY(Particle.CONNECT),2)*dt;   % particle velocity y

% Update particle's position
Particle.x(:,1)             = Particle.x(:,1) + sum(Particle.S.*VelocityX(Particle.CONNECT),2)*dt;              % particle position x
Particle.x(:,2)             = Particle.x(:,2) + sum(Particle.S.*VelocityY(Particle.CONNECT),2)*dt;              % particle position x

% Update particle's velocity gradient
Particle.Gradvelocity(:,1)  = sum(Particle.dSx.*VelocityX(Particle.CONNECT),2);                                 % particle velocity gradient xx
Particle.Gradvelocity(:,2)  = sum(Particle.dSy.*VelocityX(Particle.CONNECT),2);                                 % particle velocity gradient xy
Particle.Gradvelocity(:,3)  = sum(Particle.dSx.*VelocityY(Particle.CONNECT),2);                                 % particle velocity gradient yx
Particle.Gradvelocity(:,4)  = sum(Particle.dSy.*VelocityY(Particle.CONNECT),2);                                 % particle velocity gradient yy

% Gradvelocity = [1 2]  = [xx xy] = [vxdNx vxdNy]
%                [3 4]  = [yx yy] = [vydNx vydNy]

% Update particle's deformation gradient
% Fxx = (1+dvxxdt)*Fxx + dvxydt*Fyx
% Fxy = (1+dvxxdt)*Fxy + dvxydt*Fyy
% Fyx = dvyxdt*Fxx + (1+dvyydt)*Fyx
% Fyy = dvyxdt*Fxy + (1+dvyydt)*Fyy
Particle.defgrad_Old = Particle.defgrad;

% particle deformation gradient xx
Particle.defgrad (:,1)  = (1+Particle.Gradvelocity(:,1)*dt).*Particle.defgrad_Old(:,1) + (Particle.Gradvelocity(:,2)*dt).*Particle.defgrad_Old(:,3);
% particle deformation gradient xy
Particle.defgrad (:,2)  = (1+Particle.Gradvelocity(:,1)*dt).*Particle.defgrad_Old(:,2) + (Particle.Gradvelocity(:,2)*dt).*Particle.defgrad_Old(:,4);
% particle deformation gradient yx
Particle.defgrad (:,3)  = (Particle.Gradvelocity(:,3)*dt).*Particle.defgrad_Old(:,1) + (1+Particle.Gradvelocity(:,4)*dt).*Particle.defgrad_Old(:,3);
% particle deformation gradient yy
Particle.defgrad (:,4)  = (Particle.Gradvelocity(:,3)*dt).*Particle.defgrad_Old(:,2) + (1+Particle.Gradvelocity(:,4)*dt).*Particle.defgrad_Old(:,4);

% Update particle's volume
Particle.Jacobian = Particle.defgrad (:,1).*Particle.defgrad (:,4) - Particle.defgrad (:,2).*Particle.defgrad (:,3);
Particle.volume = Particle.Jacobian .* Particle.volume_ini;