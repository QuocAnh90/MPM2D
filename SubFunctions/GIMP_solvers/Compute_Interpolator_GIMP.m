function [Particle,Cell] = Compute_Interpolator_GIMP(Particle,Cell,Node)

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
 xnx        = xn(:,1);                                                              % nodal coordinate in X direction
 Dx         = (repmat(xp(:,1),1,Particle.Node) - xnx(Particle.CONNECT));            % distance x from point to node
 [Sx,dSx]   = GIMPshape(Dx,le(1),Particle.size(1)/2);

% 1D shape function in Y direction
 xny        = xn(:,2);                                                              % nodal coordinate in Y direction
 Dy         = (repmat(Particle.x(:,2),1,Particle.Node) - xny(Particle.CONNECT));    % distance y from point to node
 [Sy,dSy]   = GIMPshape(Dy,le(2),Particle.size(2)/2);

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
 
 function [N,dN]=GIMPshape(dX,h,lp)
%% COMPUTE BASIS FUNCTIONS
lp = 2*lp                                                                 ;% length of mp domain
c1 = ( abs(dX)< (  0.5*lp)                        )                       ;% logical array 1
c2 = ((abs(dX)>=(  0.5*lp)) & (abs(dX)<(h-0.5*lp)))                       ;% logical array 2
c3 = ((abs(dX)>=(h-0.5*lp)) & (abs(dX)<(h+0.5*lp)))                       ;% logical array 3
% BASIS FUNCTION
N1 = 1-((4*dX.^2+lp.^2)./(4*h.*lp))                                       ;% basis function according to c1
N2 = 1-(abs(dX)./h)                                                       ;% basis function according to c2
N3 = ((h+0.5*lp-abs(dX)).^2)./(2*h.*lp)                                   ;% basis function according to c3
N  = c1.*N1+c2.*N2+c3.*N3                                                 ;% basis function
% BASIS FUNCTION GRADIENT
dN1= -((8*dX)./(4*h.*lp))                                                 ;% gradient basis function according to c1
dN2= sign(dX).*(-1/h)                                                     ;% gradient basis function according to c2
dN3=-sign(dX).*(h+0.5*lp-abs(dX))./(h*lp)                                 ;% gradient basis function according to c3
dN = c1.*dN1+c2.*dN2+c3.*dN3                                              ;% gradient basis function
%-------------------------------------------------------------------------%
 