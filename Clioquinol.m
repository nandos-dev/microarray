%% Mircoarray Analysis: Effects on Clioquional on Yeast 
% Template by [Fernando Ramirez Thinh Nguyen]
% Code by [Fernando Ramirez Thinh Nguyen]
% adapted from fetch.m 
% adapted from https://www.mathworks.com/help/bioinfo/ug/working-with-geo-series-data.html
% background -- Clioquional - family of durge hydroxyquinolines, inhibit
% particular enzymes related to DNA replication. Drugs dound to have
% activity against both virla and protozoal infections 

% 00.) Analayze the microarray dataset made available by the following study 
% https://www.ncbi.nlm.nih.gov/pubmed/21504115 
% clioquinol.yeast.Li2010.pdf
%%
% 01.) download the follwing GSE1757_series_matrix.txt file from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE17nnn/GSE17257/matrix/
% index. This part is done manually. Extract the the compress txt.gz file
% using winRAR, move file to the working file directory. 

%EDA - early data analysis (exploration stage)  

gseData=bmes_downloadandparsegse_thinh_fernando('GSE17257');
get(gseData.Data); %understanding the size of the Rownames, ColNames
d = gseData.Data; % Exploring GSE data, row names and column names 
%% 

gpl_platform = gseData.Header.Series.platform_id; %saving pointer to gpl_platform  
gpl = bmes_downloadandparsegpl_thinh_fernando(gpl_platform); %obtaining metadata for the gsl_platform
gpl.ColumnNames; %outputing to function, ensuring the cell array and metadata is read

%Exploring the probsets to gene symbols, from the gplData, string comparison to the ID and Gene Symbol 
gplProbesetIDs = gpl.Data(:, strcmp(gpl.ColumnNames, 'ID')); 
geneSymbols = gpl.Data(:, strcmp(gpl.ColumnNames, 'Gene Symbol')); 
gseprobes = d.rownames;
%row_change_geneSymbol = rownames(gse.Data.Data, ':', geneSymbols); 
%the above code can be optimizing computationally by using a regex and
%saved to a variable, then variable is called. 

%mapping the GSE to GPL values, intializing a zero non-vale matrix 
MAP_GSE_GPL = zeros(numel(gseprobes),1);

%mapping of the geneSymbols same {} double as the gplProbesetIDs 
%For each gseprobe, we need to search gplprobes and use the corresponding
%gene. Doing string comparison for each of them will be too slow. Let's
%use a Map container to speed this up.

map = containers.Map(gplProbesetIDs,1:numel(gplProbesetIDs));
for i = 1:numel(gseprobes)
    if map.isKey(gseprobes{i}); MAP_GSE_GPL(i)= map(gseprobes{i});
    end
end

gsegenes = gseprobes; %make a copy, so entries not found will keep the probe name.
%genenames = gseprobes;
gsegenes(find(MAP_GSE_GPL)) = gplProbesetIDs(MAP_GSE_GPL(find(MAP_GSE_GPL))); 
%genenames(find(MAP_GSE_GPL)) = geneSymbols(MAP_GSE_GPL(find(MAP_GSE_GPL)));
d = d.rownames(':',gsegenes); 
%datamatrix = d.rownames(':',genenames);
%% Data Analysis: Determine sample groups we'll work with
% 0.3) We are often interested in comparing groups of samples. We need to look
% at the header information and decide which information for samples we can
% use to group, Usually the Header.Samples structure usually contains what we need.

%taking a dive into the header information 
samplegroups  = gseData.Header.Samples.characteristics_ch1(2,:);
unique_samplegroups = unique(samplegroups)';
%in total there are 6 samplegroups, however two uniques ones as follows 
% 1% DMSO (dimethly sulfoxide)
% 80uM CQ (tumor development in chemotherapeutic agens, anticancer drug Chloroquine)

% create logical vectirs for sample groups of interest 
Idmso = strcmpi(samplegroups,'media supplement: 1% DMSO');
Icq = strcmpi(samplegroups,'media supplement: 80 µM CQ');
%from a column indexing with logical vector 

% create a numerical vector to assign each sample to a group 1-2. As `I`
% indexed groups. 
Igroups=zeros(1,numel(samplegroups)); %initialize samplegroup size of 6 as groups  
Igroups(Idmso) = 1;
Igroups(Icq) = 2;
groupnames={'DMSO' 'CQ'};

colnames=d.colnames;
for i=1:2; colnames(Igroups==i) = groupnames(i); end
d=d.colnames(':',colnames); %this really means: "d.colnames=colnames;"
%%
%04.) Show a hiearchical clustering of samples. (Just a hiearchical
%clustering (ie a dendrogram) of samples, not a heatmap of expression
%values. 

% One idea is to only keep the genes that vary most across samples
% (ingoring sample groups.) This can be done using:
% lets try to to push this to computing to 90 variance level. 
I =genevarfilter(d,'Percentile',90);
d2 = d( I, :);

% Let's create a distance matrix between pairs of genes.
% pdist() gives a vector (to save space). If you want the symmetric matrix,
% just pass the result through squareform().
% argument changed to the spearman 
% spearman is produced as a vector, 

genedist = pdist(d2,'corr');
%'spearman'

%linkage group under each group. these results contain information about
%which two groups are combine at each branch. 
%Agglomerative hierarcvhical cluster tree. Using the average method. 
%average = UPGMA -- Unweighted average distance. 
tree = linkage(genedist,'average');

% visualize the tree, show only 20 nodes. Want to clear figure 
% Groups of genes will have a numerical id for labels.
 %bmes_fig geneclust; clf
dendrogram(tree,20,'Labels',d2.rownames);
h=gca;
h.XTickLabelRotation=45;
title('Hierarchical clustering of 20 nodes')
%distance = pdist(genedist);
%leafOrder = optimalleaforder(tree,distance) 
%%
%05.) Show a clustergram (heatmap,combined with clustering of samples and
%clustering of genes of expression values 

cg = clustergram(d2,'Standardize','Row');
%% 
%3.) Report the top 10 most different genes between the Clioquinol and
%control groups. 
[dpvals] = mattest(d(:,Idmso), d(:, Icq),'permute',1000);
[dpvals2] = mattest(datamatrix(:,Idmso), d(:, Icq),'permute',1000);
%performing two-sample t-test to evaluate differential expression of genes 
%from two experimental conditions or phenotypes, in this case it is DSMO
%and CQ medium conditions 
signif_dpvals = dpvals(dpvals(:,1) <= 0.01,:);
signif_d = d(dpvals(:,1) <= 0.01,:);
d_sig = d(dpvals(:,1) <= 0.01,:);
%taking fold change 
log2fc = log2(mean(d_sig(:,Idmso),2) ./ mean(d_sig(:,Icq),2));
%taking log2 scale to `compress scaling of values
%scatter(log2fc, -log10(signif_dpvals(:,1)), '.');
%xlabel('log_2(dsmo:CQ) media suppplements in yeast'), ylabel('-log_{10}(pvalue)');

negfc = 2.^log2fc;
negfc(negfc<1) = - 1./negfc(negfc<1);
%in order to compute the 10 most different genes between the Clioquinol and
%control groups signif_d is needed 
% Add the foldchange information to the dpvals object:
signif_dpvals=[signif_dpvals bioma.data.DataMatrix(negfc,'ColNames',{'negfc'})];
% Select the genes with pvalue<=0.01 and FC>=1.5.

I = signif_dpvals(:,'p-values')<=0.01 & abs(signif_dpvals(:,'negfc'))>=1.5;
%I = abs(signif_dpvals(:,'negfc'))>=1.5;    
%logical vector return 
dsigfc = signif_dpvals(I,:);
dsigfc = dsigfc.sortrows('p-values');
fprintf('Found %d genes with pvalue<=0.01 and FC>=1.5. Showing top 10:\n',size(dsigfc,1));
disp(dsigfc(1:10,:))

%% WRITING TO EXCEL to be used with DAVID 
I=find(signif_dpvals(:,1)<=0.01);
nsig=numel(I);
xlsdata = cell(nsig, 3); %each row will contain genesymbol,pvalue,negfc
for i=1:nsig
	gene=signif_dpvals.rownames{I(i)};
	p=signif_dpvals.double(I(i), 1);
	nfc=signif_dpvals.double(I(i), 2);
	xlsdata(i,:) = {gene p nfc};
end

xlsdata=[ {'genesymbol' 'pvalue' 'negfc'}; xlsdata]; %add the header row.
xlswrite('cq.xlsx',xlsdata,'siggenesDMSO_cq');
%%
%4.)Report the functional annotations (Go Biological Processses and KEGG
%Pathways) that are significantly different between the two groups. 

% We have found the significantly different genes between two groups. But
% what do these genes do? Are there significant differences in biological
% functions between two groups? To answer these questions, we'll make use
% of the Gene Ontology terms, which annotate each gene to one or more
% Biological Processes, Cellular Components, and Molecular Functions.

%imgo processes 
figure(1)
imshow(imread('gobp1.png'))
figure(2)
imshow(imread('gobp2.png'))
figure(3)
imshow(imread('gobp3.png'))
%KEGG pathways 
figure(4)
imshow(imread('KEGG.png'))


%% 
%5.) Discuss whether your results align with the finding reported in the
% paper 
%the results in the paper discuss gene of which fold changes were
%calculated using hte ratio of signals in C1-treated samples and DMSO
%treated controls, as shown in the analysis above. The genes that were
%found to have identical fold changes as stated in the paper and in the
%analysis above were FRE3, FET3, ENB1, ZPS1,ZRT1,ZRT3,PCA1, and SMF1. As 
%discuss in the paper the genes were upregulated fold growth as analyed in
%our excel as abs(negfc). It would be nice to mapp these genes with the
%Affymetrix probeIds, as this would eliminate the need to reference the
%geneSymbols in the code above to ensure that each ProbeID is referecing
%the correct Gene in the paper. This improvement could be added in a
%revised iteration. Likewise a p-value threshold in the paper of p<0.01 was
%considered, likewise in this analysis. For the purposes of utilizing David
%bioinformatics database, probeset IDs were used initialize to map and
%report the significant functional annotations of Go Biological Processes
%and KEGG pathways,respectively. 



