clear all

% This script runs local PCA.
 
%%
addpath(genpath('Generate_Data'),genpath('localPCA'),genpath('Sample_Results')) % add folders to path

load('compare_results.mat','data_mean','data_std','mName','d_act','D_emb','methods') % load data, i.e. results from compare.m

nManifold = length(mName); % number of manifolds tested
nMethods = length(methods); % number of methods used
nSamp = 5000;

%%
% put mean/std of other algorithms on plot for local PCA
for i = 1:nManifold
    disp(i)
    if i == 3 % generate hyperplane padded with 0s
        X = gen_plane(d_act(i),D_emb(i),nSamp,true);
    elseif i == 10 || i == 11 % generate hyperplane same dim/linearly transformed
        X = gen_plane(d_act(i),D_emb(i),nSamp,false);
    elseif i == 15 % generate hypersphere linearly transformed
        X = GenerateManifoldData(1,d_act(i)+1,nSamp);
        X = transform(X,D_emb(i));
    else % generate data (each col = 1 data pt)
        X = GenerateManifoldData(i-1,D_emb(i),nSamp); 
    end
    
    T = local_pca(X,min(d_act(i)+5,D_emb(i))); % run local PCA
    figure % plot T vs. d
    plot(0:D_emb(i),T); xlabel('d'); ylabel('Average Projection Error')
    str = strcat(num2str(d_act(i)), {'D '}, mName(i), {' in '}, num2str(D_emb(i)), 'D');
    title(str)
    hold on
    for j = 1:nMethods % add estimate from other algorithms onto plot
        errorbar(data_mean(i,j,end),0.1*j*max(T),data_std(i,j,end),data_std(i,j,end),'horizontal','*')
    end
    ylim([0 max(T)])
    plot([d_act(i) d_act(i)],[0 max(T)],'--k')
    legend(['Local PCA' methods 'Actual ID'])
    hold off
end

%%
% calculate neighborhood sizes to test
K = zeros(1,nSamp); i = 1; 
while i == 1 || K(i-1) ~= 1
    K(i) = round(nSamp/2^(i-1));
    i = i + 1;
end
K(K==0) = []; K(end) = [];
K(mod(K,2) == 1) = K(mod(K,2) == 1) + 1;

Kk = (K/2); Kk(Kk == 1) = [];

%%
% run local PCA results using different neighborhood sizes 
nM = 1:15; nM([3, 10, 15]) = []; % manifolds to test
for i = nM
    disp(i)
    if i == 3 % generate hyperplane padded with 0s
        X = gen_plane(d_act(i),D_emb(i),nSamp,true);
    elseif i == 10 || i == 11 % generate hyperplane same dim/linearly transformed
        X = gen_plane(d_act(i),D_emb(i),nSamp,false);
    elseif i == 15 % generate hypersphere linearly transformed
        X = GenerateManifoldData(1,d_act(i)+1,nSamp);
        X = transform(X,D_emb(i));
    else % generate data (each col = 1 data pt)
        X = GenerateManifoldData(i-1,D_emb(i),nSamp); 
    end
    
    T = []; Tk = [];
    for j = 1:length(K)
        disp(j)
        T(j,:) = local_pca(X,K(j),nSamp); % run local PCA
        if j <= length(Kk), Tk(j,:) = local_pca_2k(X,Kk(j),nSamp); end % run 2k local PCA
    end
    
    save(sprintf('T%d.mat',i), 'T'); % save result from local PCA
    save(sprintf('Tk%d.mat',i), 'Tk'); % save result from 2k local PCa
end

%%
% plot local PCA results using different neighborhood sizes (plot results
% from previous section)
nM = 1:15; nM([3, 10, 15]) = []; % manifolds to test
for i = nM   
    load(sprintf('T%d.mat',i)) % load results 
    load(sprintf('Tk%d.mat',i))
    
    figure % new figure for each manifold type
    str = {}; strk = {}; j = 1;
    for k = 1:length(K)
        subplot(1,2,1) % plot result for local PCA
        plot(0:D_emb(i),T(k,:))
        hold on
        
        if k <= length(Kk)
            subplot(1,2,2) % plot result for 2k local PCA
            plot(0:D_emb(i),Tk(k,:))
            strk{j} = strcat('k=',num2str(Kk(k))); % for legend
        end
        
        str{j} = strcat('k=',num2str(K(k))); % for legend
        j = j+1;
        hold on
    end
    
    subplot(1,2,2) % plot labels/formatting
    ylabel('Average Fraction of Variance Not Explained')
    ylim([0 max(Tk,[],'all')])
    plot([d_act(i) d_act(i)],[0 max(Tk,[],'all')],'--k')
    strk{end+1} = 'actual ID';
    legend(strk); xlabel('d')
    title('2k Local PCA')
    hold off
    
    subplot(1,2,1) % plot labels/formatting
    ylabel('Average Fraction of Variance Not Explained')
    ylim([0 max(T,[],'all')])
    plot([d_act(i) d_act(i)],[0 max(T,[],'all')],'--k')
    str{end+1} = 'actual ID';
    legend(str); xlabel('d')
    str = strcat(num2str(d_act(i)), {'D '}, mName(i), {' in '}, num2str(D_emb(i)), 'D');
    title([str 'Local PCA'])
    hold off
end