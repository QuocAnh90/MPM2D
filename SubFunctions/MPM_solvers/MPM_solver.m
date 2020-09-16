function [Particle] = MPM_solver(SolidModel, Cell, Node, Particle, Time, Physics)

%% Calculate the value of shape function and gradient of the shape function
[Particle,Cell] = Compute_Interpolator_MPM(Particle,Cell,Node); 

 %% Mapping from particle to nodes
[Node]=Interpolate_Particle_To_Grid(Node,Particle);

%% Update momentum and acceleration
[Node]=Integration(Node,Time,Physics,Particle);

%% Update solid particle velocity and position
[Particle] = Update_Particle(Particle,Node,Time);

%% Update particle's stress
[F_sp,V_sp,s_sp,p_sp] = Update_Stress(CModel,CModel_parameter,...
    NODES,dt,cellCount,mspoints,CONNECT,nvelo_si,dN,...
    F_sp,V_spo,m_sp,s_sp,p_sp,V_sp);