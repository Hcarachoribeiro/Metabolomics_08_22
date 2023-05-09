# Metabolomics_08_22
Scripts used to perform the statistical analysis conatained in the article "Metabolomic and proteomic profiling in bipolar disorder patients revealed potential molecular signatures related to hemostasis". Metabolomics, 18(65). https://doi.org/10.1007/s11306-022-01924-5

################## ENGLISH ########################

This Readme file comprises the proposed workflow used in our article for data processing of metabolomics and proteomics data, as well as the data integration steps and logistic regression applied to metabolomics data.

Basically, the data input in metabolomics is assumed to be already pre-processed (i.e. data after filtering, normalization, log transformation and data scaling). The data used was obtained from a HP5890 gas cromatograph coupled to a HP 5970 Series quadrupole mass spectrometer (Agilent, Santa Clara, CA, USA) and pre-proessed with MS-Dial software. For more details on the analytical side, please refer to the original article (doi number above).

#### Data after pre-processing >> PCA_and_sPLSDA_MS_Dial_VIP_final.R and PCA_and_sPLSDA_Prot_VIP_final.R

Multivariate processing from proteomics and metabolomics data, PCA from metabolomics to evaluate QCs as well. VIP selection is on the scripts, but wasn't used; Main contribuitors were used instead (next step). Output from these two scritps are the R objects for metaobolomics and proteomics used in the data integration.

#### Metabolomics and proteomics data from last step >> MultiBlock_Meta_Prot_VIP_final.R

Data integration using multiblock sPLS-DA, mixOmics package. Main Contributors (loadings to maximize covariance between Y & X) were selected as variables of interest in both datasets and exported as R objects for the next steps. Make sure the samples are the same in both datasets.

#### Multiblock-selected varibales (Proteomics and metabolomics) >> Factor_Analysis_final.R

Factor analysis was perfomed using the results from multiblock analysis. In this step the demographic and clinical data is also required to evaluate the explained variance of each factor in the dataset.

#### Multiblock-selected varibales (Proteomics and metabolomics) >> Partial2FirstOrder_Corr_Prot_Mets_Clin_final.R

Partial correlation analysis using multiblock sPLS-DA data, demographic and clinical data also needed to build correlation network. Output is the partial correlation network between proteins, metabolites and clinical/demographic data.

#### Multiblock-selected varibales (Proteomics and metabolomics) >> Bean_Plots_significants_final.R

Beanplots for both datasets. Change the input dataset(proteins/metabolites) accordingly in the script to evaluate each dataset. Output is the beanplots from each omic dataset.

#### Metabolomics R object and selected variables from multiblock analysis >> Logistic_regression_Metabolites_final.R

Logistic regression for the metabolomics dataset. Note that we used the extended metabolomics dataset here, since we analyzed more samples than the proteomics dataset. This script makes use of the variables selected in the beanplots (Snom.sign.pval.ttest) so follow the proposed order here.
Please note that in cv.glmnet function, changing the alpha values between 0 and 1 determine the type of logistic regression (ridge = 0/ LASSO = 1)
Output is the ROC curve containing the selected panel with the best potential for differentiation between the two groups.


################## PORTUGUÊS ########################

Este arquivo Readme compreende o fluxo de trabalho proposto utilizado em nosso artigo para processamento de dados de metabolômica e proteômica, bem como os passos de integração de dados e regressão logística aplicada aos dados de metabolômica.

Basicamente, a entrada de dados em metabolômica é assumida como já pré-processada (ou seja, dados após filtragem, normalização, transformação logarítmica e escalamento de dados). Os dados utilizados foram obtidos a partir de um cromatógrafo gasoso HP5890 acoplado a um espectrômetro de massa quadrupolo HP 5970 Series (Agilent, Santa Clara, CA, EUA) e pré-processados com o software MS-Dial. Para mais detalhes sobre a parte analítica,por favor consulte o artigo original (número de doi acima).

#### Dados após pré-processamento >> PCA_and_sPLSDA_MS_Dial_VIP_final.R and PCA_and_sPLSDA_Prot_VIP_final.R

Processamento multivariado de dados de proteômica e metabolômica, PCA de metabolômica para avaliar QCs também. A seleção dos VIPs está nos scripts, mas não foi utilizada; foram utilizados os principais contribuidores (próximo passo). O output desses dois scripts são os objetos R para metabolômica e proteômica usados na integração de dados.

#### Dados metabolomicos e proteomicos do passo anterior >> MultiBlock_Meta_Prot_VIP_final.R

Integração de dados usando sPLS-DA multibloco, pacote mixOmics. Os principais contribuidores (loadings para maximizar a covariância entre Y e X) foram selecionados como variáveis de interesse em ambos os conjuntos de dados e exportados como objetos R para os próximos passos. Certifique-se de que as amostras sejam as mesmas em ambos os conjuntos de dados.

#### Variáveis selecionadas na análise multibloco (Proteomica and metabolomica) >> Factor_Analysis_final.R

A análise de fatores foi realizada usando os resultados da análise multibloco. Nesta etapa, os dados demográficos e clínicos também são necessários para avaliar a variância explicada de cada fator no conjunto de dados.

#### Variáveis selecionadas na análise multibloco (Proteomica and metabolomica) >> Partial2FirstOrder_Corr_Prot_Mets_Clin_final.R

Análise de correlação parcial utilizando dados sPLS-DA multibloco, dados demográficos e clínicos também são necessários para construir uma rede de correlação. O output é a rede de correlação parcial entre proteínas, metabólitos e dados clínicos/demográficos.

#### Variáveis selecionadas na análise multibloco (Proteomica and metabolomica) >> Bean_Plots_significants_final.R

Gráficos de beanplots para ambos os conjuntos de dados. Altere o conjunto de dados de entrada (proteínas/metabólitos) de acordo no script para avaliar cada conjunto de dados. O output são os gráficos de beanplots de cada conjunto de dados omicos.

#### Conjunto de dados de metabolômica e variáveis selecionadas na análise multibloco >> Logistic_regression_Metabolites_final.R

Regressão logística para o conjunto de dados de metabolômica. Note que usamos o conjunto de dados de metabolômica estendido aqui, já que analisamos mais amostras do que o conjunto de dados de proteômica. Este script faz uso das variáveis selecionadas nos gráficos de beanplots (Snom.sign.pval.ttest), portanto siga a ordem proposta aqui.
Observe que, na função cv.glmnet, alterar os valores alpha entre 0 e 1 determina o tipo de regressão logística (ridge = 0 / LASSO = 1).
A saída é a curva ROC contendo o painel selecionado com o melhor potencial de diferenciação entre os dois grupos.


