%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <cgmeehan@alumni.caltech.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef plt_PointDensity < handle
    properties(Constant) 
        DEFAULT_BAND_WIDTH=.014;
        DEFAULT_GRID_SIZE=256;
        N_MIN=5000;
    end
    
    properties(SetAccess=private)
        deltas;
        onScale;
        xyData;
        xgrid;
        ygrid;
        wmat;
        fmat;
        fmatVector;
        ye;
        xm;
        ym;
        h;
        M;
        N;
        D;
        mins;
        maxs;
        dataBinIdxs;
    end
    
    methods
        function this=plt_PointDensity(xyData, mins, maxs, gridSize, bandWidth)
            assert(nargin>0);
            [this.N, this.D]=size(xyData);
            assert(this.D==2); %2D only
            if nargin<5
                bandWidth=this.DEFAULT_BAND_WIDTH;
                if nargin<4
                    gridSize=this.DEFAULT_GRID_SIZE;
                    if nargin<3
                        maxs=[];
                        if nargin<2
                            mins=[];
                        end
                    end
                end
            end
            needsScaling=~isempty(maxs) || ~isempty(mins);
            if isempty(maxs)
                maxs=max(xyData);
            end
            if isempty(mins)
                mins=min(xyData);
            end
            if needsScaling
                this.onScale=MatBasics.FindOnScale(xyData, mins, maxs);
                cntEdge=size(xyData, 1)-sum(this.onScale);
                if cntEdge>0
                    xyData=xyData(this.onScale, :);
                    this.N=this.N-cntEdge;
                end
                
            end
            this.mins=mins;
            this.maxs=maxs;
            this.M =gridSize;
            this.deltas=1/(this.M-1)*(maxs-mins);              
            this.xyData=xyData;
            this.h=zeros(1, this.D);
            for i =1:this.D
                if this.N<plt_PointDensity.N_MIN
                    this.h(i)=(this.N/plt_PointDensity.N_MIN)^(-1/6)*1.7...
                        *bandWidth*(maxs(i)-mins(i));
                else
                    this.h(i)=bandWidth*(maxs(i)-mins(i));
                end
            end
            this.computeWeight;
            this.computeDensity;
        end
                
        function drawJetColors(this, ax, numberOfJetColors)
            if nargin<3
                numberOfJetColors=128;
            end
            this.drawColors(ax, numberOfJetColors);
        end
        
        function H=drawContours(this, ax, percent, color, lineWidth)
            if nargin<5
                lineWidth=0.5;
                if nargin<4
                    color=[.5 .5 .6];
                    if nargin<3
                        percent=10;
                    end
                end
            end
            numLevels=floor(100/percent);
            levels=this.computeLevels(numLevels);
            [~,H]=contour(ax, this.xm, this.ym, this.fmat, levels, ...
                'k', 'color', color, 'LineStyle', '-', 'LineWidth', lineWidth);
        end
        
        function drawColors(this, ax, numberOfJetColors, subsetOfData, ...
            colorRangeStart, colorRangeEnd)
            wasHeld=ishold(ax);
            if ~wasHeld
                hold(ax);
            end
            if nargin<5
                isJetColor=true;
                if nargin<4
                    subsetOfData=[];
                    if  nargin<3
                        numberOfJetColors=64;                        
                    end
                end
            else
                isJetColor=false;
            end
            if isempty(subsetOfData)
                data=this.xyData;
            else
                if size(subsetOfData, 1)==length(this.onScale)
                    data=subsetOfData(this.onScale, :);
                else
                    onScale_=MatBasics.FindOnScale(subsetOfData, ...
                        this.mins, this.maxs);
                    data=subsetOfData(onScale_, :);
                end
            end
            z=reshape(1:this.M^2,this.M,this.M);
            eb=interp2(this.xgrid, this.ygrid, z',...
                data(:,1),data(:,2),'nearest');  %this associates each data point with its nearest grid point
            [x1,~,x3]=unique(eb);
            if isJetColor
                colors=jet(numberOfJetColors);
            else
                if nargin<6
                    colors=colorRangeStart;
                else
                    a1=mean(colorRangeEnd);
                    a2=mean(colorRangeStart);
                    if a1<a2
                        m1=max(colorRangeEnd);
                        if m1>.75
                            f=.75/m1;
                            colorRangeEnd=colorRangeEnd*f;
                        end
                    end
                    nColors=32;
                    colors=zeros(nColors,3);
                    colors(1,:)=colorRangeStart;
                    colors(nColors,:)=colorRangeEnd;
                    gap=zeros(1,3);
                    for i=1:3
                        gap(i)=colorRangeEnd(i)-colorRangeStart(i);
                    end
                    for i=2:nColors-1
                        for j=1:3
                            colors(i,j)=colors(1,j)+(i/nColors*gap(1,j));
                        end
                    end
                end
            end
            nColors=length(colors);
            try
                levels = this.computeLevels(nColors);
            catch ex
            end
            if size(data,1)<10
                color=colors(1, :);
                plot(ax, data(:,1), data(:,2), 'd',...
                    'markersize', 5, 'MarkerEdgeColor',...
                    color, 'LineStyle', 'none');
                if ~wasHeld
                    hold(ax);
                end
                return;
            end
            try
                colormap(colors);
            catch ex
                disp('huh');
            end
            densities=this.fmatVector(x1);
            lookup = plt_PointDensity.bsearch(levels,densities);
            eventColors=lookup(x3);
            
% %             usedColors=unique(eventColors);
% %             N2=length(usedColors);
% %             sz=size(data,1);
% %             if sz<10000
% %                 marker='d';
% %                 ms=2;
% %             else
% %                 marker='.';
% %                 ms=2;
% %             end
% %             for i=1:N2
% %                 colorIdx=usedColors(i);
% %                 li=eventColors==colorIdx;
% %                 plot(ax, data(li,1), data(li,2), marker,...
% %                     'markersize', ms, 'MarkerEdgeColor',...
% %                     colors(colorIdx, :), 'LineStyle', 'none');
% %             end

                % Changed this to use scatter instead of plots in loops -CS 11/15/19
                scatter(ax,data(:,1), data(:,2),5,eventColors,'filled')

            if ~wasHeld
                hold(ax);
            end
        end
    end
    
    methods(Access=private)        
        function computeWeight(this)
            this.ye=zeros(2, this.M);
            pointLL=zeros(this.N,2);  %this will be the "lower left" gridpoint to each data point
            for ii = 1:2
                this.ye(ii,:) = linspace(this.mins(ii), this.maxs(ii), this.M);
                pointLL(:,ii)=floor((this.xyData(:,ii)-this.mins(ii))./this.deltas(ii)) + 1;
            end
            pointLL(pointLL==this.M)=this.M-1;  %this avoids going over grid boundary
            %% assign each data point to its closest grid point
            [this.xgrid, this.ygrid]=meshgrid(this.ye(1,:),this.ye(2,:));
            z=reshape(1:this.M^2, this.M, this.M);
            this.dataBinIdxs=interp2(this.xgrid, this.ygrid,z',...
                this.xyData(:,1),this.xyData(:,2),'nearest');  %this associates each data point with its nearest grid point
            
            %% compute w
            Deltmat=repmat(this.deltas, this.N,1);
            shape=this.M*ones(1,2);
            wmat_=zeros(this.M, this.M);
            for ii=0:1  %number of neighboring gridpoints in 2 dimensions
                for j=0:1
                    pointm=pointLL+repmat([j ii],this.N,1);  %indices of ith neighboring gridpoints
                    pointy=zeros(this.N,2);
                    for k=1:2
                        pointy(:,k)=this.ye(k,pointm(:,k));  %y-values of ith neighboring gridpoints
                    end
                    W=prod(1-(abs(this.xyData-pointy)./Deltmat),2);  %contribution to w from ith neighboring gridpoint from each datapoint
                    wmat_=wmat_+accumarray(pointm,W,shape);  %sums contributions for ith gridpoint over data points and adds to wmat
                end
            end
            this.wmat=wmat_;
            
        end
        
        function computeDensity(this)
            Z_=zeros(1, this.D);
            Zin=cell(1, this.D);
            for i =1:this.D
                Z_(i)=min(floor(4*this.h(i)/this.deltas(i)), this.M-1);
                Zin{i}=-Z_(i):Z_(i);
            end
            phi = @(x) 1/sqrt(2*pi)*exp(-x.^2./2);
            [L_{1},L_{2}]=meshgrid(Zin{1},Zin{2});
            Phix=phi(L_{1}*this.deltas(1)./this.h(1))./this.h(1);
            Phiy=phi(L_{2}*this.deltas(2)./this.h(2))./this.h(2);
            Phimat = (Phix.*Phiy)';   
            fMat = 1/this.N*conv2(this.wmat,Phimat,'same');
            this.fmatVector=reshape(fMat,[1, this.M^2]);
            this.fmat=fMat';
            [this.xm, this.ym]=meshgrid(this.ye(1,:),this.ye(2,:));
        end
        
        function levels=computeLevels(this, numberOfLevels)
            T=sort(reshape(this.fmat, 1, this.M^this.D));
            CT=cumsum(T);
            NT=CT/CT(end);
            levels=zeros(1, numberOfLevels);
            for level=1:numberOfLevels
                idx = plt_PointDensity.bsearch(NT, level/numberOfLevels);
                levels(level)=T(idx);
            end
        end
    end
    
    methods(Static)
        function Draw(ax, data, doContours, doJetColors, reset)
            if nargin<5
                reset=true;
                if nargin<4
                    doJetColors=true;
                    if nargin<3
                        doContours=true;
                    end
                end
            end
            pb=plt_PointDensity(data);
            if reset
                cla(ax, 'reset');
            end
            try
                wasHeld=ishold(ax);
                if ~wasHeld
                    hold(ax, 'on');
                end
                if doJetColors
                    pb.drawJetColors(ax);
                end
                if doContours
                    pb.drawContours(ax);
                end
                if ~wasHeld
                    hold(ax, 'off');
                end
            catch ex
                ex.getReport
            end
        end
        
        function PlotDensity3D(ax, data, nBins, display, xLabel, yLabel, zLabel)
            if nargin<4
                display='plot';
                if nargin<3
                    nBins=64;
                end
            end
            if isempty(nBins)
                nBins=64;
            end
            if isempty(display)
                display='plot';
            end
            [R,C]=size(data);
            [D, weight, I] = plt_PointDensity.Get3D(data,nBins);
            if isequal(display, 'plot')
                jets=jet;
                cla(ax, 'reset');
                hold(ax, 'on');
                
                nClrs=size(jets,1);
                maxD=max(D(D(:)>0));
                minD=min(D(D(:)>0));
                rangeD=maxD-minD;
                d=zeros(1, R);
                for j=1:R
                    d(j)=D(I(j,1), I(j,2), I(j,3));
                end
                ratios=(d-minD)./rangeD;
                denominator=25;
                for j=1:denominator
                    ratio=j/denominator;
                    l=ratios<ratio & ratios>=(j-1)/denominator;
                    clrIdx=floor(ratio*nClrs);
                    clr2=jets(clrIdx,:);
                    plot3(ax, data(l,1), data(l,2), data(l,3), '.', ...
                        'markerSize', 1, 'lineStyle', 'none', ...
                        'markerEdgeColor', clr2, ...
                        'markerFaceColor', clr2);
                end
                view(ax, 3);
                grid(ax, 'on');
                if nargin>4
                    xlabel(ax, xLabel);
                    if nargin>5
                        ylabel(ax, yLabel)
                        if nargin>6
                            if C>2
                                zlabel(ax, zLabel);
                            end
                        end
                    end
                end
                set(ax, 'plotboxaspectratio', [1 1 1])
        
            elseif isequal(display, 'iso')
                cla(ax, 'reset');
                hold(ax, 'on');
                colormap(ax, jet);
                patch(isocaps(D,.5),...
                    'FaceColor','interp','EdgeColor','none');
                p1 = patch(isosurface(D,.5),...
                    'FaceColor',[0 .20 1],'EdgeColor','none');
                isonormals(D,p1);
                view(ax, 3);
                axis(ax, 'vis3d');
                axis(ax, 'tight')
                camlight headlight;
                lighting(ax, 'gouraud');
                grid(ax, 'on');
                set(ax, 'plotboxaspectratio', [1 1 1])
            else
                sliceomatic(weight);
            end
        end
        
        function [density, weight, idxs] = Get3D(data, nBins, bandWidth)
            if nargin<2
                nBins=64;
            end
            if nargin<3
                bandWidth=floor(nBins/5);
                if mod(bandWidth, 2)==0
                    bandWidth=bandWidth+1;
                end
            end
            x=data(:,1);
            y=data(:,2);
            z=data(:,3);
            N=numel(x);
            xBins=linspace(min(x),max(x), nBins);
            yBins=linspace(min(y),max(y), nBins);
            zBins=linspace(min(z),max(z), nBins);
            weight=zeros(nBins, nBins, nBins);
            if nargout>1
                idxs=zeros(N,3);
                for i=1:N
                    xi=find((x(i)>=xBins), 1, 'last');
                    yi=find((y(i)>=yBins), 1, 'last');
                    zi=find((z(i)>=zBins), 1, 'last');
                    weight(yi, xi, zi)=weight(yi, xi, zi)+1;
                    idxs(i,:)=[yi, xi, zi];                    
                end
            else
                for i=1:N
                    xi=find((x(i)>=xBins), 1, 'last');
                    yi=find((y(i)>=yBins), 1, 'last');
                    zi=find((z(i)>=zBins), 1, 'last');
                    weight(yi, xi, zi)=weight(yi, xi, zi)+1;
                end
            end
            density=smooth3(weight, 'gaussian', bandWidth, .8);
        end
        
        function index = bsearch(x,var)
            xLen = length(x);
            xRow= size(x,1);
            if x(1) > x(xLen)	% means x is in descending order
                if xRow==1
                    x = fliplr(x);  
                else
                    x = flipud(x);
                end
                flipped = 1;
            elseif x(1) < x(xLen)	% means x is in ascending order
                flipped = 0;
            else
                'badly formatted data. Type ''help bsearch\'')';
                return;
            end
            N=length(var);
            index=zeros(1,N);
            for i = 1:N
                low = 1;
                high = xLen;
                if var(i) <= x(low)
                    index(i) = low;
                    continue;
                elseif var(i) >= x(high)
                    index(i) = high;
                    continue;
                end
                flag = 0;
                while (low <= high)
                    mid = round((low + high)/2);
                    if (var(i) < x(mid))
                        high = mid;
                    elseif (var(i) > x(mid))
                        low = mid;
                    else
                        index(i) = mid;
                        flag = 1;
                        break;
                    end
                    if (low == high - 1)
                        break
                    end
                end
                if (flag == 1)
                    continue;
                end
                if (low == high)
                    index(i) = low;
                elseif ((x(low) - var(i))^2 > (x(high) - var(i))^2)
                    index(i) = high;
                else
                    index(i) = low;
                end
            end

            if flipped
                index = xLen - index + 1;
            end
        end
    end
end