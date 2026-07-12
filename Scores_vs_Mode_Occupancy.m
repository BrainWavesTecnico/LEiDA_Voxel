function Scores_vs_Mode_Occupancy(P,Scores_Table,Key_Modes_KC,results_dir,save_name)
% Scores_vs_Mode_Occupancy correlates the occupancy of a set of modes with a
% set of clinical/cognitive scores (partial correlation, controlling for age).
%
% P is taken directly as an input (e.g. P_original or P_harmonized from
% Save_Occupancies_Harmonize) rather than loaded from a LEiDA stats file, so
% this function does not depend on LEiDA_stats_Voxel_FracOccup_ComBat having
% been run. This makes it usable for studies with no discrete conditions to
% compare, only continuous scores to correlate with mode occupancy.
%
% INPUT:
%   P            - Fractional occupancy matrix (N_scans x length(rangeK) x rangeK(end)),
%                  e.g. P_original or P_harmonized from Save_Occupancies_Harmonize.
%   Scores_Table : .mat file with the Scores_ADNI table.
%   Key_Modes_KC : Nx2+ matrix with one row per mode to analyze, [ki c ...].
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
% Author: Joana Cabral, Tecnico, University of Lisbon, joanabcabral@tecnico.ulisboa.pt

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

figure('Color','w')

for Mode=1:N_Modes
    for Score=1:length(Selected_scores)

        if Mode==1

        disp(Scores_ADNI.Properties.VariableNames{Selected_scores(Score)})
        end

        col = Scores_ADNI{:, Selected_scores(Score)};
        valid_values = find(~isnan(col));
        n_valid_per_score(Score)=numel(valid_values);

        %[Rho(Mode,Score,1), Rho(Mode,Score,2)]=corr(col(valid_values),P_Mode(Mode,valid_values)','Type','Pearson');


        [Rho(Mode,Score,1), Rho(Mode,Score,2)] = partialcorr(col(valid_values), P_Mode(Mode,valid_values)', Age(valid_values), 'rows', 'complete');


    end

    subplot_tight(1,N_Modes+1,1+Mode,0.05)

    barh(Rho(Mode,:,1),'FaceColor',[0.84 .84 .84],'EdgeColor','none','BarWidth',.3)
    title(['Mode ' num2str(Mode)])
    hold on 
    barh(find(Rho(Mode,:,2)<0.05/N_Modes),Rho(Mode,(Rho(Mode,:,2)<0.05/N_Modes),1),'FaceColor',[0.6 .6 .6],'EdgeColor','none','BarWidth',.3)
    barh(find(Rho(Mode,:,2)<0.05/N_Modes/length(Selected_scores)),Rho(Mode,(Rho(Mode,:,2)<0.05/N_Modes/length(Selected_scores)),1),'FaceColor',Mode_colors(Mode,:),'EdgeColor','none','BarWidth',.3)
    xlim([-0.25 0.25])

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
%
