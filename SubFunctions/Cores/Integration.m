function [Node]=Integration(Node,Time,Physics,Particle)

%% Input
Nn = Node.Count;                % Total number of nodes
dt = Time.timestep;             % Time step

%% SOLVE EXPLICIT MOMENTUM BALANCE EQUATION
% COMPUTE GLOBAL NODAL ACCELERATION AND VELOCITY
Node.acceleration       = Node.inForce./Node.mass + repmat(Physics.gravity,Nn,1);
Node.velocity           = Node.momentum./Node.mass;

ID                      = Node.mass==0;
Node.acceleration(ID,:) = 0;
Node.velocity(ID,:)     = 0;

Node.velocity           = Node.velocity + dt * Node.acceleration;

%% Testing
% for n = 1:Nn
%     if Node.mass(n)>0
%         Node.velocity(n,:)      = Node.momentum(n,:)/Node.mass(n);
%         Node.acceleration(n,:)  = Node.inForce(n,:)/Node.mass(n) + [0 -10];
%         
%         Node.velocity(n,:)      = Node.velocity(n,:) + Node.acceleration(n,:)*dt;
%     end
% end

%% BOUNDARY CONDITIONS: FIX DIRICHLET BOUNDARY CONDITIONS
Node.acceleration(Node.BoundaryX,1)     =0;
Node.acceleration(Node.BoundaryY,2)     =0;
Node.velocity(Node.BoundaryX,1)         =0;
Node.velocity(Node.BoundaryY,2)         =0;


