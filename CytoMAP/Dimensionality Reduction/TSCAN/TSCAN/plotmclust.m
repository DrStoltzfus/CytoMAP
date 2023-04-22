function plotmclust(mclustobj, x, y, MSTorder, cell_size) 
    if nargin < 5
        cell_size = 6;
        if nargin < 4
            MSTorder = nan;
            if nargin < 3
                y = 2;
                if nargin < 2
                    x = 1;
                end
            end
        end
    end
    figure;
    states = size(mclustobj.clucenter,1);
    sign = ['O','v','square','<','>','+','o','x','s','d'];
    color = ['r','y','b','g','c','p','m','k','w'];
    leg = cell(1, states);
    for cid = 1:states
        temp = mclustobj.pcareduceres(mclustobj.clusterid == cid, :);
        plot(temp(:,x), temp(:, y), sign(cid),'MarkerSize', cell_size, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color(cid));
        hold on;
        leg(cid) = {['State ',num2str(cid)]};
    end
    plot(mclustobj.clucenter(:,x), mclustobj.clucenter(:, y),'Color','b','LineWidth',3);
    hold on;
    scatter(mclustobj.clucenter(:,x), mclustobj.clucenter(:, y), 55,'b','filled');
    leg(states + 1) = {'Tree'}; 
    leg(states + 2) = {'Center'};
    legend(leg);
end
