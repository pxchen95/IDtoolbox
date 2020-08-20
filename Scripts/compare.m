clear all

% This script compares the performace of various ID estimation algorithms
% on toy manifold problems across different sample sizes. 

%%
% Information on what we are testing
addpath(genpath('PCA'),genpath('MLE'),genpath('kNN'),genpath('Hein'),genpath('Generate_Data')...
    ,genpath('DANCo'),genpath('2NN')) % add relevant folders to path

mName = {'Sinusoid','Sphere','Hyperplane Padded with 0s','Strange Figure','Manifold','Helix',...
    'Manifold','Swiss Roll','Manifold','Hyperplane','Hyperplane Linearly Embedded','Moebius Band',...
    'Multivariate Gaussian','Curve','Sphere Linearly Embedded'}; % manifold shapes

d_act = [1 2 10 4 4 2 6 2 12 10 10 2 10 1 10]; % actual ID
D_emb = [3 3 100 6 8 3 36 3 72 10 100 3 10 10 100]; % extrinsic/embedding dimension
nManifold = length(d_act); % number of manifold types

nSamp = [100 1000 10000]; % sample sizes
nnSamp = length(nSamp); % number of sample sizes

nTrials = 10; % number of trials

methods = {'PCA','kNN','CD','Hein','MLE','DANCo','2NN'}; % method names
nMethods = length(methods); % number of methods

data = zeros(nManifold,nMethods,nnSamp,nTrials); % 4D matrix of ID estimates (Manifold x Method x Sample Size x Trial #)

%%
% Generate input manifold data 
X = zeros(nManifold,nTrials,max(D_emb),max(nSamp)); % manifold data (Manifold x Trial x D x Sample #)
% X_std = zeros(nManifold,nTrials); % standard deviation of manifold data (use to scale noise?)

for i = 1:nManifold % generate manifold data and calculate std of manifold data
    for j = 1:nTrials
        if i == 3 % generate hyperplane padded with 0s
            Y = gen_plane(d_act(i),D_emb(i),nSamp(end),true);
            X(i,j,1:D_emb(i),:) = Y; 
%             X_std(i,j) = min(std(Y(1:d_act(i),:),0,2));
        elseif i == 10 || i == 11 % generate hyperplane same dim/linearly transformed
            Y = gen_plane(d_act(i),D_emb(i),nSamp(end),false);
            X(i,j,1:D_emb(i),:) = Y;
%             X_std(i,j) = min(std(Y,0,2));
        elseif i == 15 % generate hypersphere linearly transformed
            Y = GenerateManifoldData(1,d_act(i)+1,nSamp(end));
            Y = transform(Y,D_emb(i));
            X(i,j,1:D_emb(i),:) = Y;
%             X_std(i,j) = min(std(Y,0,2));
        else % generate data (each col = 1 data pt)
            Y = GenerateManifoldData(i-1,D_emb(i),nSamp(end)); 
            X(i,j,1:D_emb(i),:) = Y;
%             X_std(i,j) = min(std(Y,0,2));
        end
    end
end

%%
% Run estimation algorithms
h = waitbar(0,'Please wait...'); % make wait bar
l = 1; % counter for waitbar

for i = 1:nManifold % for each manifold type
    for j = 1:nnSamp % for each sample size   
        for k = 1:nTrials % for each trial
            Y = reshape(X(i,k,1:D_emb(i),1:nSamp(j)),D_emb(i),nSamp(j));
            
            data(i,1,j,k) = dim_PCA(Y,0.01); % basic global PCA

            data(i,2,j,k) = nearneighbor(Y,4,0.01,8); % basic k nearest neighbor

            Z = GetDim(Y); 
            data(i,3,j,k) = Z(2); % correlation dimension
            data(i,4,j,k) = Z(1); % Hein (smoothed version of CD)

            data(i,5,j,k) = mledim(Y,6,20); % maximum likelihood estimator

            data(i,6,j,k) = DANCo(Y); % DANCo

            data(i,7,j,k) = twoNN(Y); % 2 nearest neighbors

            l = l + 1; % update counter for waitbar
            waitbar(l/(nManifold*nnSamp*nTrials)) % update waitbar
        end
    end    
end

close(h) % end waitbar

%%
% Plot results as a scatter plot of sample size vs. ID estimate for each manifold type
data_mean = mean(real(data),4); % mean across trials (manifold x method x sample size)
data_std = std(real(data),0,4); % standard deviation across trials

for i = 1:nManifold % make scatter plot of sample size vs ID estimate for each manifold type
    figure
    for j = 1:nMethods
        errorbar(log(nSamp),reshape(data_mean(i,j,:),1,nnSamp),reshape(data_std(i,j,:),1,nnSamp),'-o')
        hold on
    end
    plot([log(nSamp(1))-1,log(nSamp(end))+1],[d_act(i),d_act(i)],'--k')
    str = strcat(num2str(d_act(i)), {'D '}, mName(i), {' in '}, num2str(D_emb(i)), 'D w/o noise');
    title(str); xlabel('log(Sample Size)'); ylabel ('ID Estimate')
    xlim([log(nSamp(1))-1,log(nSamp(end))+1])
    legend(methods); hold off
end

%%
% Plot results as bar plot, each row = sample size, each col = grouped by
% manifold is the difference b/t true ID and est from each algorithm tested
data_diff = zeros(size(data)); % difference between estimate and true ID
for i = 1:nManifold
    d = ones(size(data(i,:,:,:)))*d_act(i);
    data_diff(i,:,:,:) = data(i,:,:,:) - d;
end

diff_mean = mean(real(data_diff),4); % mean diff across trials (manifold x method x sample size)
diff_std = std(real(data_diff),0,4); % std of diff across trials

for i = 1:nManifold % manifold names for x-axis labels
    txt = strcat(num2str(d_act(i)), {'D '}, mName(i), {' in '}, num2str(D_emb(i)), 'D');
    c = uicontrol('Style','text');
    wrap_txt = textwrap(c,txt);
    cat_txt = strjoin(wrap_txt,'\\newline');
    x_str{i} = convertCharsToStrings(cat_txt);
end

figure % make bar plot of results
j = 1; 
for i = 1:2:nnSamp
    subplot(3,1,j)
    y = abs(diff_mean(:,:,i)); % mean difference
    bar(y)
    set(gca,'xticklabel',[])
    ylabel([num2str(nSamp(i)) ' Samples'])
    hold on
    
    % ERROR BAR STUFF (1 std from mean)
    % Finding the number of groups and the number of bars in each group
    ngroups = size(y, 1);
    nbars = size(y, 2);
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    for k = 1:nbars
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*k-1) * groupwidth / (2*nbars);
        errorbar(x, y(:,k), diff_std(:,k,i), 'k', 'linestyle', 'none');
    end
    % END ERROR BAR STUFF
    
    ylim([0, max(diff_mean+diff_std,[],'all')]) %ylim([min(diff_mean,[],'all'), max(diff_mean,[],'all')])
    hold off
    j = j+1;
end
set(gca,'xticklabel',x_str)
subplot(3,1,2); legend(methods)
subplot(3,1,1); title('Absolute Difference b/t Estimate and Actual ID')