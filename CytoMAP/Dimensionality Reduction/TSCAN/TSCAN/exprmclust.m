function lpsmclust = exprmclust(data, Dreduc_opts, reduce, cluster)
    if nargin < 4
        cluster = 0;
        if nargin < 3
            reduce = true;
        end
    end
    
    if reduce == true
        [~,~,latent] = pca(normalize(data'));
        sdev = sqrt(latent');
        x = 1: min(20,size(sdev,2));
        if size(sdev,2) > 20
            sdev = sdev(1:20);
        end
        point = zeros(1,9);
        for i = 2:10
            x2 = max(0, x-i);
            [~,~,r,~,~] = regress(sdev',[ones(size(sdev));x;x2]');
            point(i-1) = sum(r.^2);
        end
        [~,optpoint] = min(point);
        pcadim = optpoint + 1;
        tmpdata = zscore(data')';
        tmppc = pca(normalize(tmpdata'));

        pcareduceres = tmpdata' * tmppc(:,1:pcadim);
        rotation = repmat([-1, 1],[1,ceil(size(pcareduceres, 2)/2)]);
        rotation = rotation(1:size(pcareduceres, 2));
        pcareduceres = pcareduceres.*rotation;
    else 
        pcareduceres = data';
    end
    
    if (cluster == 0)
        [clusterid, clunum] = Mclust(pcareduceres, Dreduc_opts);
        
    else
        clunum = length(unique(cluster));
        clusterid = cluster;
    end
    
    clucenter = zeros(clunum, size(pcareduceres,2));
    for cid = 1:clunum
        clucenter(cid, :)= mean(pcareduceres(clusterid == cid, :), 1);
    end
    
    dp = squareform(pdist(clucenter'));
    gp = graph(dp);
    dp_mst = minspantree(gp);
    
    lpsmclust.pcareduceres = pcareduceres;
    lpsmclust.MSTtree = dp_mst;
    lpsmclust.clusterid = clusterid;
    lpsmclust.clucenter = clucenter;
end
