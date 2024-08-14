clear;clc;

stored_graph;
output;

dim_per_sample = 16;
num_opt = 16+4;

plot_nodes = false;

clf;
hold on
for i = 1:(size(Points,1))
    d = ReachableVertices(((i-1)*dim_per_sample+1):i*dim_per_sample,:);
    pt = Points(i,:);
    center(i,:) = pt;
    v{i} = d;

    % d = Sol(((i-1)*num_opt+1):i*num_opt,:);
    % p(i,:) = d(1:dim_per_sample)'*v{i};
    % viol = norm(d(dim_per_sample+1:end));
    % 
    % if viol < 0.25
    %     if plot_nodes
    %         patch(v{i}(:,1),v{i}(:,2),'r','facealpha',0.001);
    %     end
    %     color(i,:) = [1 0 0];
    % else
    %     if plot_nodes
    %         patch(v{i}(:,1),v{i}(:,2),'b','facealpha',0.025);
    %     end
    %     color(i,:) = [0 1 0];
    % end
    color(i,:) = [0 0 1];
end
if plot_nodes
    scatter(center(:,1),center(:,2),100,color,'filled');
    % scatter(p(:,1),p(:,2),10,color,'filled');
end
% 
% if plot_nodes
%     x_edges = [];
%     y_edges = [];
%     for i = 1:size(Edges,1)
%         x_edges = [x_edges [center(Edges(i,1)+1,1),center(Edges(i,2)+1,1)] NaN];
%         y_edges = [y_edges [center(Edges(i,1)+1,2),center(Edges(i,2)+1,2)] NaN];
%     end
%     plot(x_edges, y_edges,'k');
% end

O_p = [];
path_plot = [];
mpc_plot = [];
start_v = [];
end_v = [];
axis([-4 4 -4 4])
axis square
mpc_N = 50;
x_ind = 1:mpc_N*4;
u_ind = (mpc_N*4+1):(mpc_N*4+(mpc_N-1)*2);

while(1)
for pt = 1:size(Path,2)
    tic
    delete(O_p)
    delete(path_plot);
    delete(mpc_plot);
    delete(start_v);
    delete(end_v);
    P = Path{pt};
    path_x = [center(P(1),1)];
    path_y = [center(P(1),2)];
    for i = 1:size(P,1)-1
        path_x = [path_x center(P(i+1),1)];
        path_y = [path_y center(P(i+1),2)];
    end
    path_plot = plot(path_x, path_y,'r','linewidth',1);
    start_v = scatter(center(P(1),1),center(P(1),2),100,'b','filled');
    end_v = scatter(center(P(end),1),center(P(end),2),100,'b','filled');
    num_traj = 10;
    mag = 1.0;
    for obs = 1:length(Obstacle_A)
    nom = lcon2vert(Obstacle_A{obs}(:,1:2), Obstacle_b{obs});
    inds = convhull(nom);
    nom = nom(inds,:)';
    Obstacle = [nom(1,:) + Obs{pt}(obs,1); nom(2,:) + Obs{pt}(obs,2)];
    O_p(obs) = patch(Obstacle(1,:),Obstacle(2,:),'r','facealpha',0.1);
    end
    

    x = MPC{pt}(x_ind);
    u = MPC{pt}(u_ind);
    x = reshape(x, 4, [])';
    u = reshape(u, 2, [])';
    
    mpc_plot = plot(x(:,1),x(:,2),'bo-','linewidth',5);

    drawnow
    val = toc;
    % pause(0.2 - val);
end
end