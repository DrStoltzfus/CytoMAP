function [label, k] = Mclust(data, Dreduc_opts)
    addpath('EmGm/EmGm');
%     [label, model, llh] = mixGaussEm(data', 3);
    
%     gmfit = fitgmdist(data,3);
%     label = cluster(gmfit,data);

% % %     Dreduc_opts = cell(4,2);
% % %     % Options for Mclust
% % %     Dreduc_opts(:,1) = {'Column',...
% % %         'CovarianceType', ...
% % %         'SharedCovariance',...
% % %         'StartCondition'};
% % %     % Defaults
% % %     Dreduc_opts(:,2) = {'2',...
% % %         'diagonal', ...
% % %         'false',...
% % %         '1'};
% % %     
% % %     dlg_title = 'Input Channel Arithmatics';
% % %     num_lines = 1;
% % %     vAnswer = inputdlg(Dreduc_opts(:,1),dlg_title,num_lines,Dreduc_opts(:,2));
% % %     
% % %     if isempty(vAnswer)
% % %         return
% % %     end
% % %     
% % %     Dreduc_opts(:,2) = vAnswer;
                
    [n, ~]= size(data);
    k = str2double(Dreduc_opts{1, 2});
    if strcmpi(Dreduc_opts{3, 2}, 'false')
        SharedCovariance = false;
    else
        SharedCovariance = true;
    end
    
    initialCond1 = randsample(1:k,n,true); 
    initialCond2 = randsample(1:k,n,true); 
    initialCond3 = 'plus'; 
    options = statset('MaxIter',1000);
    cluster0 = {initialCond1; initialCond2; initialCond3};
    
    gmfit = fitgmdist(data, k, 'CovarianceType', Dreduc_opts{2, 2},...
        'SharedCovariance', SharedCovariance, 'Start',...
        cluster0{str2double(Dreduc_opts{4, 2})},...
        'Options',options);
    label = cluster(gmfit,data);
%     for num = clusternum
%         [label, model, llh] = mixGaussEm(data', num);
%     end
end
