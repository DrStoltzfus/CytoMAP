function Plot_Region_Neighborhoods(app)

%%
smpls = fieldnames(app.data);
mlds = fieldnames(app.net);
radius = 100;
NNPlts = 40;

smpli = 5;
mdli = 11;

Ri = 7;


%% RSN Plot

fig = figure;
fig.Color = 'w';
cla

subplot(1,2,1)
cla
hold on
box off
axis equal

subplot(1,2,2)
cla
hold on
box off
axis equal


% Select model
mdl = mlds{mdli};
mdlVAR = [Constants.other_tag mdl];
%%
NR = app.net.(mdl).NR;

% Pull indeces for random neighborhoods from region Ri

Neigh = app.data.(smpls{smpli}).MFIRSN( ...
    app.data.(smpls{smpli}).MFIRSN.(mdlVAR)==Ri, ...
    {'X', 'Y', 'Z', 'NCells'});

NCellsU = unique(Neigh.NCells);
NCellsi = 1:round(numel(NCellsU)/NNPlts):numel(NCellsU);
NCellsi = NCellsU(NCellsi);

%%
Npos = zeros(numel(NCellsi),3);
cellsI = app.data.(smpls{smpli}).AllCells{:, {'X', 'Y', 'Z'}};

subplot(1,2,1)
plot(cellsI(:,1), cellsI(:,2), '.', 'color', [0.8,0.8,0.8])
    
xi = 0;
yi = 0;
for i = 1:numel(NCellsi)
    NCellsIND = find(Neigh.NCells == NCellsi(i));
    Npos(i, :) = Neigh{NCellsIND(1),{'X', 'Y', 'Z'}};
    
    % Plot the Center of the neighborhood
    x = xi+xi*radius*2;
    y = yi+yi*radius*2;
    
    % Find cells within neighborhood
    dist = sqrt(sum((Npos(i,:)-cellsI).^2, 2));
    IND = dist < radius;
    
    subplot(1,2,2)
    plot(x,y, 'xk')
    plot(cellsI(IND,1)-Npos(i,1)+x, cellsI(IND,2)-Npos(i,2)+y, '.')
    
    subplot(1,2,1)
    plot(cellsI(IND,1), cellsI(IND,2), '.')
    
    if rem(i,round(sqrt(numel(NCellsi))))==0
        yi = yi+1;
        xi = 0;
    else
        xi = xi+1;
    end
%     
end


end