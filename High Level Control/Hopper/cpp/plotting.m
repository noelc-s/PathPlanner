clear;clc;

stored_graph;
output;

dim_per_sample = 4;
num_opt = 8;

plot_nodes = false;

clf;
hold on
for i = 1:(size(Points,1))
    d = ReachableVertices(((i-1)*dim_per_sample+1):i*dim_per_sample,:);
    pt = Points(i,:);
    center(i,:) = pt;
    v{i} = d;

    % d = Sol(((i-1)*num_opt+1):i*num_opt,:);
    % p{i} = d(1:4)'*v{i};
    % viol = norm(d(5:end));
    % 
    % if viol < 0.25
    %     % patch(v{i}(:,1),v{i}(:,2),'r','facealpha',0.001);
    %     % scatter(p{i}(:,1),p{i}(:,2),10,'r','filled');
    %     color(i,:) = [1 0 0];
    % else
    %     if plot_nodes
    %         patch(v{i}(:,1),v{i}(:,2),'b','facealpha',0.025);
    %     end
    %     % scatter(p{i}(:,1),p{i}(:,2),10,'g','filled');
    %     color(i,:) = [0 1 0];
    % end
end
if plot_nodes
    scatter(center(:,1),center(:,2),100,color,'filled');
end

% x_edges = [];
% y_edges = [];
% for i = 1:size(Edges,1)
%     x_edges = [x_edges [center(Edges(i,1)+1,1),center(Edges(i,2)+1,1)] NaN];
%     y_edges = [y_edges [center(Edges(i,1)+1,2),center(Edges(i,2)+1,2)] NaN];
% end
% if plot_nodes
%     plot(x_edges, y_edges,'k');
% end

O_p = [];
path_plot = [];
start_v = [];
end_v = [];
axis([-4 4 -4 4])
axis square
while(1)
for p = 1:3:size(Path,2)
    delete(O_p)
    delete(path_plot);
    delete(start_v);
    delete(end_v);
    P = Path{p};
    path_x = [center(P(1)+1,1)];
    path_y = [center(P(1)+1,2)];
    for i = 1:size(P,1)-1
        path_x = [path_x center(P(i+1)+1,1)];
        path_y = [path_y center(P(i+1)+1,2)];
    end
    path_plot = plot(path_x, path_y,'r','linewidth',5);
    start_v = scatter(center(P(1)+1,1),center(P(1)+1,2),100,'g','filled');
    end_v = scatter(center(P(end)+1,1),center(P(end)+1,2),100,'g','filled');
    Obs = [[1 1 -1 -1] - 1.5*sin(2*3.14 * p / 1000.); [1 -1 -1 1] + 1.5*cos(2*3.14 * p / 1000.)];
    O_p = patch(Obs(1,:),Obs(2,:),'r','facealpha',0.1);
    drawnow
end
end