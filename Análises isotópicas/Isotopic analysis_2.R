#########################################################
######## Script Siar Vanessa V. Kuhnen ##################
#########################################################

# Siar = Stable Isotope Analysis in R
# Para entender melhor esse script ler antes Jackson et al. 2011 - Comparing isotopic niche widths among and within communities
# Assista também aos vídeos tutoriais em: 
# http://www.tcd.ie/Zoology/research/research/theoretical/Rpodcastsfiles/siar/siar_data_import.mp4
# http://www.tcd.ie/Zoology/research/research/theoretical/Rpodcastsfiles/siar/siar_basic.mp4
# http://www.tcd.ie/Zoology/research/research/theoretical/Rpodcastsfiles/siar/siar%20matrix%20plots.mp4

# Para não ter problemas no decorrer das análises baixe a versão 4.2.2 do SIAR disponível em: (https://github.com/AndrewLJackson/siar#siar-v422)

#############################################


rm(list = ls()) # remove os objetos antigos só por segurança
graphics.off() # fecha as janelas abertas só por segurança


library(siar) # chama o pacote do siar

setwd("D:/Vanessa/Arquivos VANESSA/aaPhD - Cap. DIETA/Análise dos dados/Análise isotópica/Siar") # Define o diretório



# ----------------------------------------------------
##    Inserindo os dados 
# -----------------------------------------------------

# Para não ter problemas no decorrer das análises de preferencia por importar arquivos com extensão .csv

# Os arquivos de dados devem sempre ter a mesma sequencia de informação nos dois arquivos: com primeiro nitrogenio, depois carbono
# Talvez o nome das colunas possa mudar, mas o tipo de informação de cada coluna deve respeitar a mesma . Mas pra garantir eu mantive os meus dados com o mesmo nome de coluna que os dados do script modelo)

# O arquivo de recursos deve ter 5 colunas nomeadas EXATAMENTE como "Means", "meand15N", "sdd15N", "meand13C", "sdd13C" - Na coluna Means vai o nome de cada recurso 

(sources <- read.csv("Dados source.csv",header=TRUE))
head(sources)

# O arquivo de consumidores deve ter 3 colunas nomeadas EXATAMENTE como "Group", "d15N", "d13C"

(consumers <- read.csv("Dados consumers.csv",header=TRUE))
head(consumers)
consumers <- as.matrix(consumers)


# Voce deve também definir concentration dependence e Fractionation inserindo arquivos com os respectivos valores, ou criando objetos iguais a zero:

#### Fractionation: 
# if you dont have any corrections make them zero (é só criar o objeto corrections<-0 )
(corrections <- read.csv("corrections.csv",header=TRUE))


#### Elemental concentration data:
# Se tiver dependencia nas concentrações insira um arquivo com esses dados: concs <- read.csv("concdepdemo.csv", header=TRUE), senão crie um objeto igual a zero
concs <- 0 




# ----------------------------------------------------
##    Rodando o modelo: 
# -----------------------------------------------------

# this line calls the SIAR model 
model1 <- siarmcmcdirichletv4(consumers, sources, corrections, concs)


# Agora dá pra fazer um biplot com os valores isotópicos dos indivíduos junto com os recursos: (precisa mudar a formatação da legenda que ficou em cima dos meus pontos)

# iso=c(2,1) é só um argumento para explicar que o eixo x deve ser a segunda coluna do arquivo dos dados que é a do carbono

siarplotdata(model1,iso=c(2,1)) 




# ---------------------------------------------------------------------------
##    Proportional contribution of each source in the consumer's diet.
# -----------------------------------------------------------------------------


# Essa função é particularmente interessante para checar se há algum recurso com uma distribuicao bi modal, mas esses gráficos dificilmente são apresentados nos trabalhos
# It will ask you which group you wish to plot the data for. You will then be asked whether you want the histograms plotted all together on one graph, or a seperate graph for each source.

siarhistograms(model1)



# ---------------------------------------------------------------------------
##    Estimates for a single group of consumers across all their sources
# -----------------------------------------------------------------------------

siarproportionbygroupplot(model1,grp=1) # D. aurita
siarproportionbygroupplot(model1,grp=2) # M. nudicaudatus


# This gets the 95% credible intervals, modes and means of the estimates. 
# It returns values for all estimated parameters (ie. the proportion of each source in the diet for each group of consumers). This value tells you how variable the consumers are within a group, after fitting the model. 
# Mode is the most likely solution (The solution with the highest probability)
# SD1G1 is the residual error associated with Isotope 1 for consumers in Group 1. (SD1 está associado com nitrogenio (lembrar que os dados do nitrogenio estavam na primeira coluna; SD2 está associado com carbono que eram os dados da segunda coluna)
# Se os valores de SD forem grandes é pq tem uma grande variação nos consumidores - isso não é bom (a variação pode vir dos sexos, das estações do ano, grau de especialização individual, etc)
(siarhdrs(model1))


# ---------------------------------------------------------------------------
##    Estimates for a single source across all groups
# -----------------------------------------------------------------------------

# aqui tb tem questão de ordem...!!!!!!!!

# Na função "siarproportionbysourceplot" apesar do argumento ser grp, que da a entender que é de grupo e que seria a espécie; na verdade o número do grp corresponde ao recurso que vc vai plotar

# the probability that the proportion of zostera in the diet of consumer group 1 is bigger than that in group 2 is then approximated by the proportion of samples that were bigger in group 1 than group 2
# Para saber qual a posição de cada recurso pra calcular a proporção: 
tabela <- as.data.frame(model1$output)
str(tabela) #  com certeza deve ter um comando melhor do que esse que eu inventei, pra não precisar ficar contando...

# Herbívoros
siarproportionbysourceplot(model1,grp=1) 
herb1 <-  model1$output[,1] >  model1$output[,9]
(herb.p1 <- sum(herb1)/length(herb1))

# Onívoros
siarproportionbysourceplot(model1,grp=2) 
Oni <-  model1$output[,2] >  model1$output[,10]
(Oni.p <- sum(Oni)/length(Oni))

# Predadores
siarproportionbysourceplot(model1,grp=3) 
Pred <-  model1$output[,3] >  model1$output[,11]
(Pred.p <- sum(Pred)/length(Pred))

# Frutos
siarproportionbysourceplot(model1,grp=4) 
fruit <-  model1$output[,4] >  model1$output[,12]
(fruit.p <- sum(fruit)/length(fruit))

# Vertebrados
siarproportionbysourceplot(model1,grp=5) 
vert <-  model1$output[,5] >  model1$output[,13]
(vert.p <- sum(vert)/length(vert))

# Detritívoros
siarproportionbysourceplot(model1,grp=6) 
detri <-  model1$output[,6] >  model1$output[,14]
(detri.p <- sum(detri)/length(detri))


# ---------------------------------------------------------------------------
##    Correlation between sources - matrix plots
# -----------------------------------------------------------------------------

# These matrix plots show how the estimated dietary proportions are correlated with each other. (ou seja... é interessante para dizer se para manter o mesmo padrão isotópico a espécie teria que aumentar o consumo dos dois itens, ou se tem que "escolher" aumentar o consumo entre um dos dois...
# A diagonal nada mais é do que o mesmo resultado da função que ja foi feita acima: siarhistograms(model1); Os valores abaixo da diagonal são os coeficiente de correlação; E acima da diagonal são só os scaterplots da correlação
# Um alto valor negativo de correlação significa que quando o modelo tenta aumentar a proporção de um ele acaba diminuindo a proporção do outro (isso pq no final a soma sempre tem q dar 1!!!)
# E um alto valor positivo significa que o aumento na contribuição de um dos recursos implica no aumento da contribuição do outro tb

# It will ask you which group you wish to run the analysis for. 
siarmatrixplot(model1)




# ------------------------------------------------------------
##        Population-Level Metrics (Layman metrics) 
# ------------------------------------------------------------


#  Arguments:
# x,y (Bivariate data given as vectors x and y)
# O x é o carbono que na nossa planilha está na segunda coluna, e o Y é o nitrogenio que está na nossa planilha na primeira coluna

#  Values:
# dN_range -> Assuming y is delta Nitrogen
# dC_range -> Assuming x is delta Carbon
# hull  -> Contains the area of the convex hull around the data points defined by x and y (hull$TA); the coordinates for plotting of the convex hull (hull$xcoords, hull$ycoords); and the index address of the points in x and y which define the convex hull (hull$chI)
# CD -> Mean distance to centroid 
# MNND -> Mean Nearest Neighbour Distance
# SDNND -> Standard Deviation of the Nearest Neighbour Distance

# Para entender melhor o que cada valor gerado significa, dar uma lida no artigo Jackson et al 2012: Population-Level Metrics of Trophic Structure Based on Stable Isotopes and Their Application to Invasion Ecology_PlosOne. 
# O script completo da função está em: https://github.com/AndrewLJackson/siar/blob/master/R/laymanmetrics.R



# ---  D. aurita:
# Para aurita os dados estão nas linhas 1:19
(x.a <- consumers[1:19,2]) # Lembrar q os dados de carbono entraram na segunda coluna!
(y.a <- consumers[1:19,1]) # Lembrar q os dados de nitrogenio entraram na primeira coluna!
(Layman.aurita <- laymanmetrics(x.a,y.a))


# ---  M. nudicaudatus:
# Para meta os dados estão nas linhas 20:33
(x.m <- consumers[20:33,2]) # Lembrar q os dados de carbono entraram na segunda coluna!
(y.m <- consumers[20:33,1]) # Lembrar q os dados de nitrogenio entraram na primeira coluna!
(Layman.meta <- laymanmetrics(x.m,y.m))





# ------------------------------------------------------------
##        Observações finais 
# ------------------------------------------------------------


# save your model to a file name of your choice
save(model1,file="model.rdata")


# If you want to access the raw data output from the model you can get it...
# This line just returns the first 10 rows (and all the columns) of the data because its very large. The histograms and density plots generated above are simply histograms of each column in this file which represents the posterior density draws.
model1$output[1:10,]

