clear;clc;

stored_graph;
output;

dim_per_sample = 16;
num_opt = 16+4;

result = yaml.loadFile("/home/noelcs/repos/ThinkingThroughThings/High Level Control//Hopper/cpp/config/params.yaml");

plot_edges = result.log_edges;
u_max = result.MPC.tau_max;

%%% Beizer
dt = result.MPC.dt;
gamma = 2;
order = 2*gamma-1; % minimal curve
m = 1;

H = Bezier.H(order, dt);
D = Bezier.D(gamma,order, dt);
Z = Bezier.Z(order, dt);
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
B = H_vec*inv(D)';
tau = linspace(0,dt);
A_x = [1 0; -1 0; 0 1; 0 -1];
b_x = [3;3;3;3];
Lf = 0;
Lg = 1;
K = [-1 -1];
e_bar = 0;
Delta_vec = Bezier.Delta_vec(m, order, gamma);
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
D_vec = Delta_vec*H_vec;
%%%

clf;
hold on
for i = 1:(size(Points,1))
    d = ReachableVertices(((i-1)*dim_per_sample+1):i*dim_per_sample,:);
    pt = Points(i,:);
    center(i,:) = pt;

    if plot_edges
        xbar = pt([1 3])';
        f_xbar = 0;
        g_xbar = 1;
        [F, G] = Bezier.F_G(A_x, b_x, H, m, xbar, f_xbar, g_xbar, gamma,{eye(4)},Lg,Lf,e_bar,K,u_max);
        Vert = cddmex('extreme',struct('A',[D_vec(1:2,:); F],'B',[xbar;G],'lin',1:2));
        Vert = Bezier.Poly.conv((D_vec(3:4,:)*Vert.V')');
        v{i} = Vert;
        patch(v{i}(:,1),pt([2])+v{i}(:,1)*0,v{i}(:,2),'g','facealpha',0.025);
        xbar = pt([2 4])';
        f_xbar = 0;
        g_xbar = 1;
        [F, G] = Bezier.F_G(A_x, b_x, H, m, xbar, f_xbar, g_xbar, gamma,{eye(4)},Lg,Lf,e_bar,K,u_max);
        Vert = cddmex('extreme',struct('A',[D_vec(1:2,:); F],'B',[xbar;G],'lin',1:2));
        Vert = Bezier.Poly.conv((D_vec(3:4,:)*Vert.V')');
        v{i} = Vert;
        patch(pt([1])+v{i}(:,1)*0,v{i}(:,1),v{i}(:,2),'r','facealpha',0.025);
    end
    
    % if plot_edges
    %     d = Sol(((i-1)*num_opt+1):i*num_opt,:);
    %     p(i,:) = d(1:dim_per_sample)'*v{i};
    %     viol = norm(d(dim_per_sample+1:end));
    % 
    %     if viol < 0.25
    %         if plot_edges
    %             patch(v{i}(:,1),v{i}(:,2),'r','facealpha',0.001);
    %         end
    %         color(i,:) = [1 0 0];
    %     else
    %         if plot_edges
    %             patch(v{i}(:,1),v{i}(:,2),'b','facealpha',0.025);
    %         end
    %         color(i,:) = [0 1 0];
    %     end
    % else
        color(i,:) = [0 0 1];
    % end
end
if plot_edges
    scatter3(center(:,1),center(:,2),center(:,3),100,repmat([0 1 0],size(center,1),1),'filled');
    scatter3(center(:,1),center(:,2),center(:,4),100,repmat([1 0 0],size(center,1),1),'filled');
    % scatter(p(:,1),p(:,2),10,color,'filled');
end

if plot_edges
    x_edges = [];
    y_edges = [];
    for i = 1:size(Edges,1)
        x_edges = [x_edges [center(Edges(i,1)+1,1),center(Edges(i,2)+1,1)] NaN];
        y_edges = [y_edges [center(Edges(i,1)+1,2),center(Edges(i,2)+1,2)] NaN];
    end
    plot(x_edges, y_edges,'k');
end

O_p = [];
path_plot = [];
mpc_plot = [];
start_v = [];
end_v = [];
axis([-4 4 -4 4])
axis square
mpc_N = result.MPC.N;
x_ind = 1:mpc_N*4;
u_ind = (mpc_N*4+1):(mpc_N*4+(mpc_N-1)*2);

% while(1)
    for pt = 1:size(Path,2)
        tic
        delete(O_p)
        delete(path_plot);
        delete(mpc_plot);
        delete(start_v);
        delete(end_v);
        P = Path{pt}+1;
        path_x = [center(P(1),1)];
        path_y = [center(P(1),2)];
        for i = 1:size(P,1)-1
            path_x = [path_x center(P(i+1),1)];
            path_y = [path_y center(P(i+1),2)];
        end
        path_plot = plot(path_x, path_y,'r','linewidth',1);
        start_v = scatter(center(P(1),1),center(P(1),2),100,'c','filled');
        end_v = scatter(center(P(end),1),center(P(end),2),100,'c','filled');
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

        Bezier_x = [];
        Bezier_y = [];
        for i = 1:size(x,1)-1
            Xi_x = B*[x(i,[1 3])'; x(i+1,[1 3])'];
            Xi_y = B*[x(i,[2 4])'; x(i+1,[2 4])'];
            Bezier_x = [Bezier_x reshape(Xi_x,2,[])*Z(tau)];
            Bezier_y = [Bezier_y reshape(Xi_y,2,[])*Z(tau)];
        end

        mpc_plot(1) = plot(x(:,1),x(:,2),'bo','linewidth',5);
        mpc_plot(2) = plot(Bezier_x(1,:),Bezier_y(1,:),'b','linewidth',5);

        drawnow
        val = toc;
        pause(0.2  - val);
    end
% end