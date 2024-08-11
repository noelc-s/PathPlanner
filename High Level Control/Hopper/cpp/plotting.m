clear;clc;

output;

dim_per_sample = 5;
num_opt = 8;

Obs = [1 1 -1 -1; 1 -1 -1 1];

clf;
hold on
for i = 1:(size(O,1)/dim_per_sample)
    i
    d = O(((i-1)*dim_per_sample+1):i*dim_per_sample,:);
    center(i,:) = d(1,:);
    v{i} = d(2:end,:);

    d = V(((i-1)*num_opt+1):i*num_opt,:);
    p{i} = d(1:4)'*v{i};
    viol = norm(d(5:end));
    
    if viol < 0.25
        patch(v{i}(:,1),v{i}(:,2),'r','facealpha',0.001);
        % scatter(p{i}(:,1),p{i}(:,2),10,'r','filled');
        color(i,:) = [1 0 0];
    else
        patch(v{i}(:,1),v{i}(:,2),'b','facealpha',0.025);
        % scatter(p{i}(:,1),p{i}(:,2),10,'g','filled');
        color(i,:) = [0 1 0];
    end
end
scatter(center(:,1),center(:,2),100,color,'filled');

x_edges = [];
y_edges = [];
for i = 1:size(E,1)
    x_edges = [x_edges [center(E(i,1)+1,1),center(E(i,2)+1,1)] NaN];
    y_edges = [y_edges [center(E(i,1)+1,2),center(E(i,2)+1,2)] NaN];
end
plot(x_edges, y_edges,'k');

for i = 1:size(P,1)-1
    plot([center(P(i)+1,1),center(P(i+1)+1,1)],[center(P(i)+1,2),center(P(i+1)+1,2)],'r','linewidth',5)
end

patch(Obs(1,:),Obs(2,:),'r','facealpha',0.1);
axis equal