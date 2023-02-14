function Figure=Plot_Final(Particle,Node,Cell)

Count = Particle.Count;
x_sp = Particle.x;
r1_sp = Particle.r1;
r2_sp = Particle.r2;
v_ssp = Particle.velocity;
s_sp = Particle.stress;
LOC = Node.x;
d_sp = Particle.displacement;

% x_corner1              = zeros(Count,2);
% x_corner2              = zeros(Count,2);
% x_corner3              = zeros(Count,2);
% x_corner4              = zeros(Count,2);
 
displacement = zeros(Count,1);
% velocity = zeros(Count,1);
% strain = zeros(Count,1);
for sp=1:Count
%     displacement(sp) = sqrt(d_sp(sp,1)^2+d_sp(sp,2)^2);
%     velocity(sp) = sqrt(v_ssp(sp,1)^2+v_ssp(sp,2)^2);
%     velocity(sp) = v_ssp(sp,1);
%     velocity(sp) = v_ssp(sp,2);
    stress(sp) = s_sp(1,sp);
end

%  for sp=1:Count
%  x_corner1(sp,:) = x_sp(sp,:) - r1_sp(sp,:) - r2_sp(sp,:);      % Position of corner 1
%  x_corner2(sp,:) = x_sp(sp,:) + r1_sp(sp,:) - r2_sp(sp,:);
%  x_corner3(sp,:) = x_sp(sp,:) + r1_sp(sp,:) + r2_sp(sp,:);
%  x_corner4(sp,:) = x_sp(sp,:) - r1_sp(sp,:) + r2_sp(sp,:);
%  end

Figure=figure
% set(Figure, 'visible','off');
sz = 10;
color = stress;
% color = b_sp(:,1);
scatter(x_sp(:,1),x_sp(:,2),sz,color,'filled');
hold on
% for sp=1:Count
%     x_cor = [x_corner1(sp,1) x_corner2(sp,1) x_corner3(sp,1) x_corner4(sp,1) x_corner1(sp,1)];
%     y_cor = [x_corner1(sp,2) x_corner2(sp,2) x_corner3(sp,2) x_corner4(sp,2) x_corner1(sp,2)];
%     plot(x_cor,y_cor,'b')
% end
grid on
axis([0,max(LOC(:,1)),0,max(LOC(:,2))]);
set(gca,'xtick',[0:10*Cell.size(1):max(LOC(:,1))]);
set(gca,'ytick',[0:10*Cell.size(2):max(LOC(:,2))]);
h=colorbar;
colormap(jet(256))
% set(h, 'ytick', [0:0.2:1.2]);
%     caxis([0 0.0001]);