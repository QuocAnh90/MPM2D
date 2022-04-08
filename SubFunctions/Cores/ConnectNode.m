function [Cell,Particle] = ConnectNode(Cell, Node, Particle,Version)

%% Switch version
switch Version
         case 'MPM'
            Particle.Node  = 4;                             % Total number of nodes
            
         case 'GIMP'
            Particle.Node  = 16;                            % Total number of nodes
end

%% Input
pCount          = Particle.Count;           % Total number of particle

%% Generate variables
Cell.CONNECT        = zeros(Cell.Count,Particle.Node);      % store nodal index of a cell
Cell.Particle       = cell(Cell.Count,1);                   % store particle index 
Particle.CONNECT    = zeros(pCount,Particle.Node);          % store nodal index
Particle.Elems      = zeros(pCount,1);                      % Cell index where store particles

Particle.S          = zeros(Particle.Node,pCount);          % value of shape function
Particle.dSx        = zeros(Particle.Node,pCount);          % value of shape function grandient in X
Particle.dSy        = zeros(Particle.Node,pCount);          % value of shape function grandient in Y

%% Calculate cell indexes

switch Version
         case 'MPM'
            for c = 1:Cell.Count
             % COmpute nodal index of element 'c'
             bottomleft = c + floor(c/(Node.CountX-1));
             bottmright = bottomleft+1; 
             topright   = bottmright+Node.CountX; 
             topleft    = bottomleft+Node.CountX;         
             Cell.CONNECT(c,:)  = [bottomleft bottmright topright topleft];       % Store nodal index
              
%                      4(x)------(x)3
%                        |        |
%                        |   o    |
%                        |        |
%                      1(x)------(x)2
            end
            
         case 'GIMP'
            for c = 1:Cell.Count
             % COmpute nodal index of element 'c'
             bottomleft = c + floor(c/(Node.CountX-1));
             bottmright = bottomleft + 1; 
             topright   = bottmright + Node.CountX; 
             topleft    = bottomleft + Node.CountX;
             
             node1      = bottomleft - Node.CountX - 1;
             node2      = node1 + 1;
             node3      = node2 + 1;
             node4      = node3 + 1;
             node5      = bottomleft - 1;
             node8      = bottmright + 1;
             node9      = topleft - 1;
             node12     = topright + 1;
             node13     = topleft + Node.CountX - 1;
             node14     = node13 + 1;
             node15     = node14 + 1;
             node16     = node15 + 1;
             
             Cell.CONNECT(c,:)  = [node1 node2 node3 node4 node5 bottomleft bottmright node8 node9 topright topleft node12 node13 node14 node15 node16];       % Store nodal index
             
%                     13(x)------(x)14-----(x)15-----(x)16
%                        |        |         |        |
%                        |        |         |        |
%                        |        |         |        |
%                      9(x)------(x)10-----(x)11-----(x)12             
%                        |        |         |        |
%                        |        |    o    |        |
%                        |        |         |        |
%                      5(x)------(x)6------(x)7------(x)8 
%                        |        |         |        |
%                        |        |         |        |
%                        |        |         |        |
%                      1(x)------(x)2------(x)3------(x)4 
            end
end


end
