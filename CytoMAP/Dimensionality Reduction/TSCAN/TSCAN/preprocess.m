function predata = preprocess(data)
    if nargin < 8
        cvcutoff = 0;
        if nargin < 7
            minexpr_percent = 0;
            if nargin < 6
                minexpr_value = 0;
                if nargin < 5
                    pseudocount = 1;
                    if nargin < 4
                        logbase = 2;
                        if nargin < 3
                            takelog = true;
                            if nargin < 2
                                clusternum = 0;
                            end
                        end
                    end
                end
            end
        end
    end
    
    if (takelog) 
        data = log(data + pseudocount) / log(logbase);
    end
    predata = data(mean(data > minexpr_value, 2) > minexpr_percent & std(data, 0, 2) ./ mean(data,2) > cvcutoff, :);
    if (clusternum ~= 0)
        clures = linkage(pdist(predata),'single');
        cluster_result = cluster(clures,'maxclust',clusternum);
        aggdata= zeros(clusternum,size(predata,2));
        for i = 1: clusternum
            aggdata(i,:) = mean(predata(cluster_result==i,:),1);
        end
        predata = aggdata;
    end
end

         