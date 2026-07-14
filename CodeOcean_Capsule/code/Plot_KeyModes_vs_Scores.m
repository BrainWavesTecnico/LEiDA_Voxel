function Plot_KeyModes_vs_Scores(P,Scores_Table,Key_Modes_KC,results_dir,save_name)
% Plot_KeyModes_vs_Scores plots and exports the correlation (partial
% correlation, controlling for age) between the occupancy of a selected set
% of key modes and a set of clinical/cognitive scores.
%
% This is the figure-generation counterpart to Scores_vs_Mode_Occupancy.m:
% that function scans the entire pyramid of modes and saves the resulting
% statistics (no figure, no dependency on which modes you'll eventually want
% to look at); this function then renders one bar plot per selected mode
% (e.g. the Key_Modes_KC returned by Choose_Relevant_Modes.m) and exports a
% CSV of the correlations.
%
% INPUT:
%   P            - Fractional occupancy matrix (N_scans x length(rangeK) x rangeK(end)),
%                  e.g. P_original or P_harmonized from Save_Occupancies_Harmonize.
%   Scores_Table : .mat file with the Scores_ADNI table.
%   Key_Modes_KC : Nx2+ matrix with one row per mode to plot, [ki c ...].
%                  ki is the POSITION of the clustering solution in rangeK
%                  (i.e. P's 2nd dimension index), NOT the literal number of
%                  clusters - e.g. if rangeK = 2:20, ki=1 means K=2 clusters.
%   results_dir  : Directory where the figure/CSV/mat outputs are saved.
%   save_name    : Output .mat filename for the correlation results.
%
% NOTE: The set of score columns used (Genetics/Biomarkers/Cognitive_functions
% indices below) is hardcoded for the ADNI Scores_ADNI table used in Campo et
% al.; adapt these indices for a different scores table.
%
% Author: Joana Cabral, University of Lisbon, joanabcabral@tecnico.ulisboa.pt

load(Scores_Table,'Scores_ADNI')

Mode_colors=[137 207 240; 227 120 91; 250 221 107; 207 225 185] ./ 256;

% Sort the centroids by the occupancy in controls
N_Modes=size(Key_Modes_KC,1);

for Mode = 1:N_Modes

    k=Key_Modes_KC(Mode,1);
    c=Key_Modes_KC(Mode,2);
    P_Mode(Mode,:) = P(:,k,c);
end

Genetics=[20];
Biomarkers=[21 23 24 25 26 27];
Cognitive_functions=[28 29  38 56 57 34 35 36 37 30 31 32 33  40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 ];

Selected_scores=flip([Genetics Biomarkers Cognitive_functions]);

Age=Scores_ADNI.AGE_AT_SCAN;

Rho=(zeros(N_Modes,length(Selected_scores),2));

% Compute correlations for every mode/score up front, so the shared x-axis
% range below can be based on the true maximum across the whole figure.
for Score=1:length(Selected_scores)
    disp(Scores_ADNI.Properties.VariableNames{Selected_scores(Score)})

    col = Scores_ADNI{:, Selected_scores(Score)};
    valid_values = find(~isnan(col));
    n_valid_per_score(Score)=numel(valid_values);

    for Mode=1:N_Modes
        [Rho(Mode,Score,1), Rho(Mode,Score,2)] = partialcorr(col(valid_values), P_Mode(Mode,valid_values)', Age(valid_values), 'rows', 'complete');
    end
end

% Shared x-axis range (same scale in every subplot) based on the largest
% absolute correlation found anywhere in the figure, with a small margin.
max_abs_rho = max(abs(Rho(:,:,1)), [], 'all');
xlim_range = [-1 1] * (max_abs_rho + 0.02);

alpha_mode   = 0.05 / N_Modes;
alpha_strict = 0.05 / N_Modes / length(Selected_scores);

figure('Color','w')

for Mode=1:N_Modes

    subplot_tight(1,N_Modes+1,1+Mode,0.05)

    % Color each bar by significance level within a SINGLE barh() call, so
    % every bar has the same width. Overlaying separate barh() calls on
    % sparse/uneven subsets of categories (as before) makes barh compute
    % 'BarWidth' relative to each call's own (uneven) category spacing,
    % which is why the highlighted bars were rendering thicker than the rest.
    bar_colors = repmat([0.84 0.84 0.84], length(Selected_scores), 1);
    sig_mode   = Rho(Mode,:,2) < alpha_mode;
    sig_strict = Rho(Mode,:,2) < alpha_strict;
    bar_colors(sig_mode,:)   = repmat([0.6 0.6 0.6], sum(sig_mode), 1);
    bar_colors(sig_strict,:) = repmat(Mode_colors(Mode,:), sum(sig_strict), 1);

    b = barh(Rho(Mode,:,1), 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', .3);
    b.CData = bar_colors;

    title(['Mode ' num2str(Mode)])
    xlim(xlim_range)

    box off

    if Mode ==1
        names  = Scores_ADNI.Properties.VariableNames(Selected_scores);
        labels = strcat(names, ' (N=', arrayfun(@num2str, n_valid_per_score, 'UniformOutput', false), ')');
        set(gca,'YTick',1:length(Selected_scores),'YTickLabel',labels,'TickLabelInterpreter','none','FontSize',8)
    else
        set(gca,'YTick',[],'FontSize',8)
    end
end

%% CREATE CSV table
% Assumes: Rho(mode, variable, 1) = r, Rho(mode, variable, 2) = p-value
% and: labels = 36x1 cell array of variable names

nVars = numel(Selected_scores);

% Build header
header = {'Variable'};
for m = 1:N_Modes
    header{end+1} = sprintf('Mode %d (Pearson r)', m);
    header{end+1} = sprintf('Mode %d (p-value)', m);
end

% Build table as cell array of strings
T = cell(nVars+1, 1 + N_Modes*2);
T(1,:) = header;

for v = 1:nVars
    T{v+1, 1} = labels{v};
    for m = 1:N_Modes
        r = Rho(m, v, 1);
        p = Rho(m, v, 2);
        T{v+1, 2*m}   = sprintf('%.3f', r);
        T{v+1, 2*m+1} = sprintf('%.2e', p);  % forces scientific notation
    end
end

% After building T (keeping header row in place)
header_row = T(1,:);
data_rows = T(2:end,:);
data_rows_flipped = flipud(data_rows);
T = [header_row; data_rows_flipped];

% Write to CSV as strings (avoids Numbers/Excel truncating small p-values)
fid = fopen(fullfile(results_dir, 'BraVe_correlations.csv'), 'w');
for row = 1:size(T,1)
    fprintf(fid, '%s', strjoin(T(row,:), ','));
    fprintf(fid, '\n');
end
fclose(fid);

save(fullfile(results_dir, save_name),'P_Mode','Age','Biomarkers','Cognitive_functions','Genetics','Scores_ADNI','Selected_scores','Rho')
