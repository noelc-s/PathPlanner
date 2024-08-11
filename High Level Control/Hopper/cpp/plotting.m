
output;

dim_per_sample = 5;
num_opt = 8;

Obs = [1 1 -1 -1; 1 -1 -1 1];

clf;
hold on
for i = 1:(size(O,1)/dim_per_sample)
    i
    d = O(((i-1)*dim_per_sample+1):i*dim_per_sample,:);
    p = d(1,:);
    v{i} = d(2:end,:);

    d = V(((i-1)*num_opt+1):i*num_opt,:);
    p = d(1:4)'*v{i};
    viol = sum(d(5:end));
    
    if viol < 0.2
        scatter(p(:,1),p(:,2),100,'r','filled');
        patch(v{i}(:,1),v{i}(:,2),'r','facealpha',0.001);
        scatter(p(:,1),p(:,2),100,'r','filled');
    else
        scatter(p(:,1),p(:,2),100,'g','filled');
        patch(v{i}(:,1),v{i}(:,2),'b','facealpha',0.1);
        scatter(p(:,1),p(:,2),100,'g','filled');
    end
end
patch(Obs(1,:),Obs(2,:),'r','facealpha',0.1);
axis equal