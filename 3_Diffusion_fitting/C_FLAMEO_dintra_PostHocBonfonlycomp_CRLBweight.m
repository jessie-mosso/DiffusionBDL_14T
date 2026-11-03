clear;clc;close all; 

%% ============================================================
%  Group-level analysis — INDEPENDENT GROUPS (between-subjects)
%  *** BONFERRONI VERSION (for comparison with FDR) ***
%  Each "week" = a different cohort of animals.
%  Within-subject fit error propagated (FLAMEO / weighted GLM style)
%% ============================================================


%% --- 1. LOAD DATA ---------------------------------------------------
% betaFile = 'parameters.xlsx';   % point estimates
% varFile  = 'variances.xlsx';    % parallel table of fit VARIANCES
%                                 % (diagonal of covariance for the
%                                 %  parameter of interest, e.g. c1 or c2)
% 
% B = readmatrix(betaFile);
% V = readmatrix(varFile);

load([pwd '/FLAMEO_input_indiv_CRLBweightedfit/dsticks.mat'])
load([pwd '/FLAMEO_input_indiv_CRLBweightedfit/dsticks_var.mat'])

dsticks(dsticks==0)=NaN; 
dsticks_var(dsticks_var==0)=NaN; 

B=dsticks;
V=dsticks_var;


metabolites = {'Gln','glu','ins','tau','tNAA','tCr'};
weeks       = [0 4 6];            % 3 INDEPENDENT GROUPS
nMetab      = numel(metabolites);
nGroup      = numel(weeks);

%% --- 2. RESHAPE INTO 3D: animals x group(week) x metabolite -------
nAnimal = size(B,1);
Bmat = nan(nAnimal, nGroup, nMetab);
Vmat = nan(nAnimal, nGroup, nMetab);
for m = 1:nMetab
    cols = (m-1)*nGroup + (1:nGroup);
    Bmat(:,:,m) = B(:, cols);
    Vmat(:,:,m) = V(:, cols);
end

%% --- 3. FLAME-1 PRECISION-WEIGHTED GROUP MEAN (per group) ---------
groupMean = nan(nGroup, nMetab);
groupSE   = nan(nGroup, nMetab);
betweenVar= nan(nGroup, nMetab);

for m = 1:nMetab
    for g = 1:nGroup
        y  = Bmat(:,g,m);
        wv = Vmat(:,g,m);
        ok = ~isnan(y) & ~isnan(wv);
        y  = y(ok);  wv = wv(ok);
        if numel(y) < 2, continue; end
        [mu, se, tau2] = flame1_mean(y, wv);
        groupMean(g,m) = mu;
        groupSE(g,m)   = se;
        betweenVar(g,m)= tau2;
    end
end

%% --- 4. OVERALL GROUP EFFECT (one-way, weighted) ------------------
omnibus = table();
for m = 1:nMetab
    [grpID, yVal, wVar] = deal([]);
    for g = 1:nGroup
        y  = Bmat(:,g,m);  wv = Vmat(:,g,m);
        ok = ~isnan(y) & ~isnan(wv);
        grpID = [grpID; repmat(g, sum(ok), 1)];
        yVal  = [yVal;  y(ok)];
        wVar  = [wVar;  wv(ok)];
    end
    if numel(unique(grpID)) < 2, continue; end

    tbl = table(yVal, categorical(grpID, 1:nGroup, ...
                strcat("wk",string(weeks))), wVar, ...
                'VariableNames',{'beta','week','wvar'});

    mdl = fitlm(tbl, 'beta ~ week', 'Weights', 1./tbl.wvar);
    avTerms = anova(mdl);
    idx = strcmp(avTerms.Properties.RowNames,'week');

    omnibus.Metabolite(m) = string(metabolites{m});
    omnibus.F(m)          = avTerms.F(idx);
    omnibus.df1(m)        = avTerms.DF(idx);
    omnibus.df2(m)        = avTerms.DF(strcmp(avTerms.Properties.RowNames,'Error'));
    omnibus.p(m)          = avTerms.pValue(idx);

    fprintf('\n=== %s : overall group effect ===\n', metabolites{m});
    disp(avTerms);
end

%% --- 5. PAIRWISE COMPARISONS BETWEEN GROUPS -----------------------
pairs = nchoosek(1:nGroup, 2);
pw = table();
row = 0;
for m = 1:nMetab
    for k = 1:size(pairs,1)
        g1 = pairs(k,1);  g2 = pairs(k,2);
        y1 = Bmat(:,g1,m); v1 = Vmat(:,g1,m);
        y2 = Bmat(:,g2,m); v2 = Vmat(:,g2,m);
        ok1 = ~isnan(y1)&~isnan(v1); ok2 = ~isnan(y2)&~isnan(v2);
        y1=y1(ok1); v1=v1(ok1); y2=y2(ok2); v2=v2(ok2);
        if numel(y1)<2 || numel(y2)<2, continue; end

        [m1,s1] = flame1_mean(y1, v1);
        [m2,s2] = flame1_mean(y2, v2);

        diff = m1 - m2;
        seD  = sqrt(s1^2 + s2^2);
        z    = diff / seD;
        df   = (s1^2+s2^2)^2 / ( s1^4/(numel(y1)-1) + s2^4/(numel(y2)-1) );
        p    = 2*tcdf(-abs(z), df);

        row = row+1;
        pw.Metabolite(row) = string(metabolites{m});
        pw.Comparison(row) = string(sprintf('wk%d - wk%d', weeks(g1), weeks(g2)));
        pw.Diff(row)       = diff;
        pw.SE(row)         = seD;
        pw.t(row)          = z;
        pw.df(row)         = df;
        pw.p(row)          = p;
    end
end

%% --- 6. MULTIPLE COMPARISONS CORRECTION (BONFERRONI) --------------
% (a) Omnibus: Bonferroni across metabolites (keep or set = p; see notes)
% omnibus.p_Bonf = min(omnibus.p * height(omnibus), 1);
omnibus.p_Bonf = omnibus.p;   % treat metabolites as independent

% (b) Pairwise: Bonferroni ONLY across the 3 pairs WITHIN each metabolite
pw.p_Bonf = nan(height(pw),1);
uMetab = unique(pw.Metabolite, 'stable');
for i = 1:numel(uMetab)
    rows = strcmp(pw.Metabolite, uMetab(i));   % the (up to 3) pairs
    nThis = sum(rows);                         % number of comparisons here
    pw.p_Bonf(rows) = min(pw.p(rows) * nThis, 1);  % multiply by n, cap at 1
end

fprintf('\n========= OMNIBUS (overall group effect) =========\n');
disp(omnibus);
fprintf('\n========= PAIRWISE COMPARISONS (Bonferroni within metabolite) =========\n');
disp(pw);

%% --- 7. PLOT: individual datapoints + group result + sig bars -----
figure('position',[100,100,2000,200])
for m = 1:nMetab
%     subplot(2, ceil(nMetab/2), m); hold on;
    subplot(1, 6, m); hold on;

    allY = [];
    for g = 1:nGroup
        y  = Bmat(:,g,m);  y = y(~isnan(y));
        allY = [allY; y];
        x  = repmat(g, numel(y),1) + 0.12*(rand(numel(y),1)-0.5);
        scatter(x, y, 25, [.6 .6 .6], 'filled', 'MarkerFaceAlpha',0.5);
    end

    errorbar(1:nGroup, groupMean(:,m), groupSE(:,m), ...
             'o','LineWidth',2,'Color',[0.2 0.4 0.8], ...
             'MarkerFaceColor',[0.2 0.4 0.8],'LineStyle','none');

    title(metabolites{m}); ylabel('Parameter');
    xlim([0.5 nGroup+0.5]); xticks(1:nGroup);
    xticklabels(strcat("wk",string(weeks)));

    % --- significance bars + stars (BONFERRONI-corrected pairwise) ---
    rowsM = find(strcmp(pw.Metabolite, metabolites{m}));
    sigRows = rowsM(pw.p_Bonf(rowsM) < 0.05);

    if ~isempty(allY)
        yTop  = max(allY);
        yBase = yTop + 0.04*range(allY);
        yStep = 0.07*range(allY);
    else
        yBase = 1; yStep = 0.1;
    end

    for ii = 1:numel(sigRows)
        r = sigRows(ii);
        cmp = char(pw.Comparison(r));
        nums = sscanf(cmp, 'wk%d - wk%d');
        g1 = find(weeks == nums(1));
        g2 = find(weeks == nums(2));

        yBar = yBase + (ii-1)*yStep;
        plot([g1 g1 g2 g2], ...
             [yBar-0.2*yStep, yBar, yBar, yBar-0.2*yStep], ...
             'k', 'LineWidth', 1.2);

        p = pw.p_Bonf(r);
        if     p < 0.001, stars = '***';
        elseif p < 0.01,  stars = '**';
        else              stars = '*';
        end
        text(mean([g1 g2]), yBar + 0.15*yStep, stars, ...
             'HorizontalAlignment','center', ...
             'FontSize', 12, 'FontWeight','bold');
    end

    if ~isempty(sigRows) && ~isempty(allY)
        ylim([min(allY)-0.05*range(allY), ...
              yBase + numel(sigRows)*yStep + 0.1*range(allY)]);
    end
end
sgtitle('Individual animals (grey) + group result (blue); * Bonferroni p<0.05');
for met=1:6
    for gr=1:3
        tableprismSE((met-1)*12+(gr-1)*3+1)=groupMean(gr,met);
        tableprismSE((met-1)*12+(gr-1)*3+2)=groupSE(gr,met);
    end 
end 

for met=1:6
    for gr=1:3
        tableprismindiv(:,(met-1)*4+1:(met-1)*4+3)=dsticks(:,(met-1)*3+1:(met-1)*3+3);
    end 
end 

tableprismSE(tableprismSE==0)=NaN;
tableprismindiv(tableprismindiv==0)=NaN;

%Gln w4-w6
(groupMean(3,1)-groupMean(2,1))./groupMean(2,1)*100
%Glu w0-w4
(groupMean(2,2)-groupMean(1,2))./groupMean(1,2)*100
%mIns w0-w6
(groupMean(3,3)-groupMean(1,3))./groupMean(1,3)*100
%Tau w0-w4
(groupMean(2,4)-groupMean(1,4))./groupMean(1,4)*100




%% ============================================================
%  LOCAL FUNCTIONS
%% ============================================================
function [mu, se, tau2] = flame1_mean(y, wv)
    tau2 = max(var(y) - mean(wv), 0);
    for it = 1:100
        w   = 1 ./ (wv + tau2);
        mu  = sum(w.*y) / sum(w);
        Q   = sum(w .* (y - mu).^2);
        num = Q - (numel(y) - 1);
        den = sum(w) - sum(w.^2)/sum(w);
        tau2new = max(num/den, 0);
        if abs(tau2new - tau2) < 1e-10, tau2 = tau2new; break; end
        tau2 = tau2new;
    end
    w  = 1 ./ (wv + tau2);
    mu = sum(w.*y) / sum(w);
    se = sqrt(1 / sum(w));
end
