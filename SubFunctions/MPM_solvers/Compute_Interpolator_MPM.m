function [Particle,Cell] = Compute_Interpolator_MPM(Particle,Cell,Node)

%% Input
xp              = Particle.x;                       % Position of particle
le              = Cell.size;                        % Size of cell
CellCount       = Cell.Count;                       % Totala number of cells
xn              = Node.x;                           % Position of nodes

%% Output
% Particle.S     % value of shape function
% Particle.dS    % Value of shape function gradient

%% Calculate the nodal index interacting with particle
 % Cell index of particle p
 Particle.Elems = floor(xp(:,1)/le(1)+1)+(Node.CountX-1)*(floor(xp(:,2)/le(2))); 
 
 % Nodal index of particle p
 Particle.CONNECT    = Cell.CONNECT(Particle.Elems,:); 
 
%% Calculate value of shape function
% 1D shape function in X direction
 xnx        = xn(:,1);   % nodal coordinate in X direction
 Dx         = (repmat(xp(:,1),1,Particle.Node) - xnx(Particle.CONNECT));  
 [Sx,dSx]   = linearshape(Dx,le(1));

% 1D shape function in Y direction
 xny        = xn(:,2);   % nodal coordinate in Y direction
 Dy         = (repmat(Particle.x(:,2),1,Particle.Node) - xny(Particle.CONNECT));  
 [Sy,dSy]   = linearshape(Dy,le(2));

% shape function in 2D
 Particle.S     = Sx.*Sy;
 Particle.dSx   = dSx.*Sy;
 Particle.dSy   = dSy.*Sx;
 
 % Build B matrix
 i            = 1:2:Particle.Node*2-1;
 j            = i+1;
 Particle.B(1,i,:) = Particle.dSx';
 Particle.B(2,j,:) = Particle.dSy';
 Particle.B(3,i,:) = Particle.dSy';
 Particle.B(3,j,:) = Particle.dSx';
 
 %% Compute Cell.Particle: index of particles in each cell (active cell)
 for c =1:CellCount
     id_p = find(Particle.Elems==c);
     Cell.Particle{c}=id_p;
 end
 
 
 function [N,dN]=linearshape(dx,Lx)
    c1  = (abs(dx)<=Lx);
    N1 = 1-abs(dx)/Lx; 
    dN1 = -sign(dx)/Lx;
    N = N1.*c1; 
    dN = c1.*dN1;
