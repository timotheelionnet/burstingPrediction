function threshInteractions = dispInteractions(M,burstingData, tfList,intList,threshMax,...
    geneSymbolToDisplay1,geneSymbolToDisplay2)

%% plot stats for cherry picked example
i1 = find(ismember(tfList,geneSymbolToDisplay1));
j1 = find(ismember(intList,geneSymbolToDisplay1));
i2 = find(ismember(tfList,geneSymbolToDisplay2));
j2 = find(ismember(intList,geneSymbolToDisplay2));

disp(['TF ',geneSymbolToDisplay1,' has ',num2str(M(i1,j2)),...
    ' evidences of interaction with interactor ',geneSymbolToDisplay2]);

disp(['TF ',geneSymbolToDisplay2,' has ',num2str(M(i2,j1)),...
    ' evidences of interaction with interactor ',geneSymbolToDisplay1]);

%% plot histograms of the number of interaction evidences
[n,x] = hist(M(:),0:1:max(M(:))); 
figure('Name','Frequency of interaction evidence'); 
plot(x,n);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
xlabel('evidence');
ylabel('Number of interactions');

%% plot histogram of number of interactions per TF

[nChar,xChar] = hist(sum(M(ismember(tfList,burstingData.TFnames),:),2),0:1:500); 
[nNonChar,xNonChar] = hist(sum(M(~ismember(tfList,burstingData.TFnames),:),2),0:1:500); 

figure('Name','Number of interactions per TF');
hold;
plot(xChar,cumsum(nChar)/sum(nChar),'DisplayName','Characterized TFs');
plot(xNonChar,cumsum(nNonChar)/sum(nNonChar),'DisplayName','NonCharacterized TFs');
xlabel('Number of interactions');
ylabel('Cumulative Fraction of TFs');

%% plot histogram of number of unique interactions per TF
[nChar,xChar] = hist(sum(M(ismember(tfList,burstingData.TFnames),:)>0,2),0:1:500); 
[nNonChar,xNonChar] = hist(sum(M(~ismember(tfList,burstingData.TFnames),:)>0,2),0:1:500); 

figure('Name','Number of unique interactions per TF');
hold;
plot(xChar,cumsum(nChar)/sum(nChar),'DisplayName','Characterized TFs');
plot(xNonChar,cumsum(nNonChar)/sum(nNonChar),'DisplayName','NonCharacterized TFs');
xlabel('Number of unique interactions');
ylabel('Cumlative Fraction of TFs');

%% prepare a stats table
summation = sum(M>0);
threshInt(:,1) = 1:threshMax;
TFintCounts(:,1) = 1:numel(tfList);

for n = 1:threshMax
    m = (find(summation>=n));
    threshInt(n,2) = numel(m);
    threshInt(n,3) = sum(sum(M(:,m)));
    
    for q = 1:numel(tfList)
        TFintCounts(q,n+1) = sum(M(q,m));
    end
end

threshInteractions = array2table(threshInt,...
    'variableNames',{'thresh', 'interactors', 'interactions'});

%% plot histogram of # of interactions / # of interactors as a function of
% threshold
figure('Name','Number of interactions/interactors as a function of threshold'); 
hold;
plot(threshInteractions.thresh,threshInteractions.interactors,'-o',...
    'DisplayName','Number of interactors');
plot(threshInteractions.thresh,threshInteractions.interactions,'-o',...
    'DisplayName','Number of interactions');
xlabel('Minimal number of interactions per interactor');
legend;














 