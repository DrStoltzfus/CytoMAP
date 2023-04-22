%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <cgmeehan@alumni.caltech.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef Supervisors < handle
    properties(Constant)
        DEV_UNIT_LIMIT=4;
        BORDER_DEV_UNIT_LIMIT=.2;
        VERBOSE=false;
        COLOR_NEW_SUBSET=[153 153 179];
    end
    
    properties
        compareByMeans=false;
        density;
    end
    properties(SetAccess=private)
        ids;
        labelMap;
        mdns;
        mads;
        means;
        stds;
        cnts;
        N;
        xLimit;
        yLimit;
        labels;
    end
    
    methods
        
        function this=Supervisors(labels, labelMap, embedding, ax)
            this.labelMap=labelMap;
            this.ids=unique(labels);
            this.labels=labels;
            N=length(this.ids);
            this.cnts=zeros(N,1);
            this.mdns=zeros(N,2);
            this.means=zeros(N,2);
            this.mads=zeros(N,2);
            this.stds=zeros(N,2);
            for i=1:N
                id=this.ids(i);
                l=labels==id;
                this.cnts(i)=sum(l);
                this.mdns(i,:)=median(embedding(l,:));
                this.mads(i,:)=mad(embedding(l,:),1);
                this.means(i,:)=mean(embedding(l,:));
                this.stds(i,:)=std(embedding(l,:),1);
            end
            this.N=N;
            if nargin<4 || isempty(ax) || ~ishandle(ax)
                mx=max(embedding);
                mn=min(embedding);
                this.xLimit=[mn(1) mx(1)];
                this.yLimit=[mn(2) mx(2)];
            else
                this.xLimit=xlim(ax);
                this.yLimit=ylim(ax);
            end
        end
        
        function [name, color, label]=getNameByMedian(this, mdn)
            [~, iMeHearty]=pdist2(this.density.clusterMdns, mdn, ...
                'euclidean', 'smallest', 1);
            name=this.density.clusterNames{iMeHearty};
            color=this.density.clusterColors{iMeHearty};
            label=this.density.clusterLabels{iMeHearty};
        end        
    end
    
    methods(Static, Access=private)
        function [density, numClusts, clustIds]=Cluster(data, mins, maxs)
            %labels=Clumps(embedding, this.mdns,this.ids, 16);
            if nargin<3
                [numClusts, clustIds, density]=Density.ClusterVeryHigh(data);
            else
                [numClusts, clustIds, density]=Density.ClusterVeryHigh(...
                    data, mins, maxs);
            end
        end
    end
    
    methods
        function matchWithClusters(this, data, density, ...
                numClusters, clusterIds)
            if nargin>3
                this.density=density;
            end          
            clusterLabels=cell(1,numClusters);
            clusterNames=cell(1,numClusters);
            clusterColors=cell(1,numClusters);
            %[~,I]=pdist2(this.mdns, embedding, 'euclidean', 'smallest', 1);
            %labels2=this.ids(I);
            if ~this.compareByMeans || isempty(this.means)
                avgs=this.mdns;
                if isprop(this, 'mads')
                    devs=this.mads;
                else
                    devs=[];
                end
            else
                avgs=this.means;
                devs=this.stds;
            end
            if any(this.cnts<20)
                avgs(this.cnts<20, 1)=this.xLimit(2)*5;
                avgs(this.cnts<20, 2)=this.yLimit(2)*5;
            end
            cluAvgs=zeros(numClusters, 2);
            cluDevs=zeros(numClusters, 2);
            for i=1:numClusters
                l=clusterIds==i;
                if Supervisors.VERBOSE
                    sum(l)
                end
                if this.compareByMeans
                    cluAvgs(i,:)=mean(data(l,:));
                    cluDevs(i,:)=std(data(l,:));
                else
                    cluAvgs(i,:)=median(data(l,:));
                    cluDevs(i,:)=mad(data(l,:), 1);
                end
            end
            hasDevs=~isempty(devs);
            [D, I]=pdist2(avgs, cluAvgs, 'euclidean', 'smallest', 1);
            labels_=zeros(size(data, 1), 1);
            reChecks={};
            for i=1:numClusters
                labelIdx=I(i);
                label=this.ids(labelIdx);
                if label==0
                    label=0-i;
                else
                    key=num2str(label);
                    clusterLabels{i}=label;
                    clusterNames{i}=this.labelMap.get(key);
                    clusterColors{i}=this.labelMap.get([key '.color']);
                end
                l=clusterIds==i;
                if Supervisors.VERBOSE
                    sum(l)
                    this.labelMap.get(num2str(label))
                end
                if hasDevs
                    devDist=MatBasics.DevDist(cluAvgs(i,:), cluDevs(i,:));
                    if any(D(i)>=devDist*Supervisors.DEV_UNIT_LIMIT)
                        reChecks{end+1}=struct('clustId', i, 'count',...
                            sum(l), 'label', label, 'labelIdx', labelIdx);
                        label=0-i;
                    end
                end
                labels_(l)=label;
            end
            if hasDevs
                N_=length(reChecks);
                while N_>0 
                    changes=[];
                    for i=1:N_
                        clustId=reChecks{i}.clustId;
                        label=reChecks{i}.label;
                        labelIdx=reChecks{i}.labelIdx;
                        closestLabelIdxs=labels_==label;
                        if Supervisors.VERBOSE
                            sum(closestLabelIdxs)
                            this.labelMap.get(num2str(label))
                            reChecks{i}
                        end
                        if any(closestLabelIdxs)
                            %Does this cluster with no label match
                            %   sit on the border of one with the 
                            %   closest label from the supervisor?
                            unlabeledClusterIdxs=clusterIds==clustId;
                            borderDistance=min(pdist2(...
                                data(unlabeledClusterIdxs, :), ...
                                data(closestLabelIdxs,:), 'euclidean', ...
                                'smallest', 1));
                            supervisorDevDistance=MatBasics.DevDist(...
                                avgs(labelIdx,:), devs(labelIdx,:));
                            limit=supervisorDevDistance*...
                                Supervisors.BORDER_DEV_UNIT_LIMIT;
                            if borderDistance<=limit
                                changes(end+1)=i;
                                labels_(unlabeledClusterIdxs)=label;
                            end
                        end
                    end
                    if isempty(changes)
                        break;
                    else
                        reChecks(changes)=[];
                    end
                    N_=length(reChecks);                    
                end
                newSubsetIds=unique(labels_(labels_<0));
                newSubsets=length(newSubsetIds);
                for i=1:newSubsets
                    clustId=0-newSubsetIds(i);
                    clusterLabels{clustId}=newSubsetIds(i);
                    clusterNames{clustId}=['New subset #' num2str(i)];
                    color_=Supervisors.COLOR_NEW_SUBSET+clustId;
                    if any(color_<0)
                        color_(color_<0)=0;
                    end
                    clusterColors{clustId}=num2str(color_);
                end
                N_=length(this.ids);
                remainder=data(labels_==0,:); 
                labels2=zeros(size(remainder,1),1);
                for i=1:N_
                    if ~any(find(labels_==this.ids(i),1))
                        label=this.ids(i);
                        if Supervisors.VERBOSE
                            this.labelMap.get(num2str(label))
                        end
                        [D2, I2]=pdist2(avgs(i,:), remainder, ...
                            'euclidean', 'smallest', 1);
                        pt=avgs(i,:)+devs(i,:);
                        devDist=pdist2(avgs(i,:), pt);
                        if any(D2<devDist*Supervisors.DEV_UNIT_LIMIT)
                            labels2(D2<devDist)=label;
                        end
                    end
                end
                if any(labels2)
                    labels_(labels_==0)=labels2;
                end
            else
                newSubsets=0;
            end
            this.density.setLabels(labels_, newSubsets, clusterNames, ...
                clusterColors, cluAvgs, cluDevs, clusterLabels);
        end
        
        function [labels, labelMap]=supervise(this, data, ...
                computeClusters, doHtml)
            if nargin<4
                doHtml=false;
            end
            if computeClusters
                [dns, numClusters, clusterIds]=Supervisors.Cluster(data);
                this.matchWithClusters(data, dns, numClusters, clusterIds);
            end              
            labels=this.density.labels;
            labelMap=java.util.Properties;
            ids_=unique(labels);
            N_=length(ids_);
            for i=1:N_
                putInMap(ids_(i));
            end
            
            function putInMap(id) 
                key=num2str(id);
                keyColor=[key '.color'];
                if id==0
                    if doHtml
                        name='<font color="#608080">background</font>';
                    else
                        name='\color[rgb]{0.36 0.36 0.5}\bf\itbackground';
                    end
                    color='92 92 128';
                elseif id<0
                    clustId=0-id;
                    nm=this.density.clusterNames{clustId};
                    if doHtml
                        name=['<font color="#9292AA">' nm ' ?</font>'];
                    else
                        name=['\color[rgb]{0.6 0.6 0.7}\bf\it' nm  ' ?'];
                    end
                    color=this.density.clusterColors(clustId);
                else
                    name=this.labelMap.get(key);
                    if doHtml
                        if String.Contains(name, '^{')
                            name=strrep(name, '^{', '<sup>');
                            name=strrep(name, '}', '</sup>');
                        end
                    end
                    color=this.labelMap.get(keyColor);
                end
                labelMap.put(key,name);
                labelMap.put(keyColor, color);
            end
            
        end
        
        function drawClusterBorders(this, ax)
            wasHeld=ishold(ax);
            if ~wasHeld
                hold(ax, 'on');
            end
            N_=length(this.density.clusterColors);
            for i=1:N_
                clr=(str2num(this.density.clusterColors{i})/256)*.85;
                gridEdge(this.density, true, i, clr, ax, .8, '.', '-', .5);
                if Supervisors.VERBOSE
                    str2num(this.density.clusterColors{i})
                    clr
                    disp('ok');
                end
            end
            if ~wasHeld
                hold(ax, 'off');
            end
        end

    end
end