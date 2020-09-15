function [Node,Cell] = Grid_Generation(Node,Cell,x_min,x_max,y_min,y_max)
%% Generate the structured grid
% Input
NN = [Node.CountX Node.CountY]; le = Cell.size;
% NN(1): number of nodes in X direction
% NN(2): number of nodes in Y direction
% le(1): element size in X direction
% le(2): element size in X direction

%% Generate the coordinates of nodes and cells
LOC                     = zeros(NN(1)*NN(2),2);             % nodal coordinate
LOCC                    = zeros((NN(1)-1)*(NN(2)-1),2);     % element centroid coordinate

LOCX                    = [0:NN(1)-1]'*le(1);               % Location of all nodes in X direction
LOCY                    = [0:NN(2)-1]'*le(2);               % Location of all nodes in Y direction

LOCCX                   = [0:(NN(1)-1)-1]'*le(1)+le(1)/2;   % Location of cells in X direction
LOCCY                   = [0:(NN(2)-1)-1]'*le(2)+le(2)/2;   % Location of cells in X direction

for i=1:NN(2)
    LOC((1+NN(1)*(i-1)):(NN(1)*(i-1)+NN(1)),1) = LOCX;      % generate the X node position in LOC
end

for i=1:NN(2)
    LOC((NN(1)*(i-1))+1:NN(1)*i,2) = LOCY(i);               % generate the Y node position in LOC
end

for i=1:NN(2)-1
    LOCC((1+(NN(1)-1)*(i-1)):((NN(1)-1)*(i-1)+(NN(1)-1)),1) = LOCCX;        % generate the X element position in LOCC
end

for i=1:NN(2)-1
    LOCC(((NN(1)-1)*(i-1))+1:(NN(1)-1)*i,2) = LOCCY(i);                     % generate the Y element position in LOCC
end

%% Detect boudary condition
% nfbcx: number of boundary nodes in X direction
% nfbcy: number of boundary nodes in Y direction
% fbcx: index of all boundary nodes in X direction
% fbcy: index of all boundary nodes in Y direction

%% Nodal variables generation
Node.Count          = NN(1)*NN(2);                  % total number of nodes
Node.x              = LOC;                          % Nodal Coordinate
Node.mass           = zeros(Node.Count,1);          % Nodal Mass
Node.momentum       = zeros(Node.Count,2);          % Nodal Momentum
Node.inForce        = zeros(Node.Count,2);          % Nodal Internal force
Node.eForce         = zeros(Node.Count,2);          % Nodal External force
Node.traction       = zeros(Node.Count,2);          % Nodal Traction
Node.acceleration   = zeros(Node.Count,2);          % Nodal Acceleration
Node.velocity       = zeros(Node.Count,2);          % Nodal Velocity

% Compute nodal index in the boundary
[nfbcx,nfbcy,fbcx,fbcy]=Compute_Boundary_Nodes(Node.Count,LOC,x_max,x_min,y_max,y_min);
Node.Count_BoundaryX =  nfbcx;          % number of boundary nodes in X direction
Node.Count_BoundaryY =  nfbcy;          % number of boundary nodes in Y direction
Node.BoundaryX = fbcx;                  % nodal index in the boundary  in X direction
Node.BoundaryY = fbcy;                  % nodal index in the boundary  in Y direction

%% Cell variables generation
Cell.Count      = (NN(1)-1)*(NN(2)-1);          % total number of elements
Cell.x          = LOCC;
