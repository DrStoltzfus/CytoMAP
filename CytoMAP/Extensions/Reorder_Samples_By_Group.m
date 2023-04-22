function Reorder_Samples_By_Group(app)
    %% This function re-orders the samples in the order defned by the group[ number metadata
    smpls = fieldnames(app.data);
    groups = zeros(1, numel(smpls));
    for smpli = 1:numel(smpls)
        if contains('MetaData', fieldnames(app.data.(smpls{smpli})))
            if contains('Group', fieldnames(app.data.(smpls{smpli}).MetaData))
                grp = app.data.(smpls{smpli}).MetaData.Group;
                if iscell(grp)
                    groups(smpli) = grp{1};
                else
                    groups(smpli) = grp;
                end
            else
                groups(smpli) = 0; 
            end
        else
            groups(smpli) = 0;
        end
    end
    [~, INDgrp]= sort(groups);
    app.data = orderfields(app.data, INDgrp);
    app.DataN.Items = fieldnames(app.data);
    app.DataN.Value = app.DataN.Items{1};

end
