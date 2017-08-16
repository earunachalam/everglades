% plot system from csv

clf

plotstyle_voronoi = 1;

for i = 499%0:49
    disp([num2str(i) ' starting'])
    
    d = csvread(['../reaction_system/u_1.' num2str(i) '.csv'],1,0);
    x = d(:,2);
    y = d(:,3);
    conc = d(:,1);
    
    if plotstyle_voronoi
        [v,c] = voronoin([x y]);
        for j = 1:length(c) 
            patch(v(c{j},1),v(c{j},2),conc(j),'EdgeColor',[0.5 0.5 0.5]);
        end
    else
        scatter(x,y,100,conc,'filled')
    end
    
    axis equal
    axis off
    xlim([min(x) max(x)])
    ylim([min(y) max(y)])
    pause(0.5);
end