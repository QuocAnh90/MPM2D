function [Cell,Particle] = ConnectNode(Cell, Node, Particle)

%% Input
pCount          = Particle.Count;           % Total number of particle

%% Generate variables
% Cell.CONNECT   = cell(Cell.Count,1);        % store nodal index
Cell.CONNECT   = zeros(Cell.Count,4);       % store nodal index
Cell.Particle  = cell(Cell.Count,1);        % store particle index 
% Particle.Node  = 4*ones(pCount,1);          % Total number of nodes interacting with a particle
% Particle.S     = cell(pCount,1);            % value of shape function
% Particle.dS    = cell(pCount,1);            % Value of shape function gradient
Particle.Node  = 4;                         % Total number of nodes
Particle.S     = zeros(pCount,4);           % value of shape function
Particle.dSx   = zeros(pCount,4);           % value of shape function
Particle.dSy   = zeros(pCount,4);           % value of shape function
Particle.Elems = zeros(pCount,1);           % Cell index where store particles
% Particle.CONNECT   = cell(pCount,1);        % store nodal index
Particle.CONNECT   = zeros(pCount,4);       % store nodal index

%% Calculate cell indexes
for c = 1:Cell.Count
         % COmpute nodal index of element 'c'
         bottomleft = c + floor(c/(Node.CountX-1));
         bottmright = bottomleft+1; 
         topright   = bottmright+Node.CountX; 
         topleft    = bottomleft+Node.CountX;         
%          Cell.CONNECT{c}  = [bottomleft bottmright topright topleft];       % Store nodal index
         Cell.CONNECT(c,:)  = [bottomleft bottmright topright topleft];       % Store nodal index
        end
end
