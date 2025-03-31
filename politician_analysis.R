
# Real data ---------------------------------------------------------------

rm(list=ls())
gc()

# Libraries ---------------------------------------------------------------

library(truncnorm)
library(mnormt)
library(psych)
library(progress)
library(invgamma)
library(mcclust)
library(ggplot2)
library(mcclust.ext)
library(reshape2)
library(viridis)
library(plot.matrix)
library(sams)
library(gtools)
library(readxl)
library(dplyr)
library(msm)
library(salso)
library(ggalluvial)

source("gibbs_cpo.R")


# Data --------------------------------------------------------------------

data <- as.data.frame(read_excel("senate_votes_may.xls"))

data[data==9|data==0] <- NA
data[data==6] <- 0

data1 <- data[,c(1:47)]
data2 <- data[,c(1:4, 48:82)]

elim1 <- elim2 <- c()

for(i in 1:nrow(data)){

  if(sum(is.na(data1[i,]))>30){

    elim1 <- c(elim1, i)
  }
  if(sum(is.na(data2[i,]))>28){

    elim2 <- c(elim2, i)

  }

}

data1 <- data1[-c(elim1), ]
data2 <- data2[-c(elim2), ]

data1[data1==7] <- NA


dd1 <- data1[,5:ncol(data1)]
rownames(dd1) <- 1:nrow(dd1)
colnames(dd1) <- 1:ncol(dd1)
dd2 <- data2[,5:ncol(data2)]
rownames(dd2) <- 1:nrow(dd2)
colnames(dd2) <- 1:ncol(dd2)
par(mar=c(5.1,4.1,4.1,4.1))
plot(as.matrix(dd1), main="Votes of the U.S. House of Representative",
     breaks=c(-0.5,0.5,1.5), xlab="Votes", ylab="Congressmen/Congresswomen")
plot(as.matrix(dd2), main="Votes of the U.S. Senate",
     breaks=c(-0.5,0.5,1.5), xlab="Votes", ylab="Senators")

# Senators -----------------------------------------------------------------

bin2 <- as.matrix(data2[,c(5:39)])
head(bin2)

data = bin2
n <- nrow(data)
p <- ncol(data)
rownames(data) <- c(1:n)
colnames(data) <- c(1:p)

plot(data, xaxt="n", yaxt="n")


sigma2 <- 0.1
tau2 <- 1.5
Z <- data*1000 + matrix(rnorm(n*p, 0 , sd=sqrt(sigma2)), ncol=p, nrow=n)

for(i in 1:n){
  for(j in 1:p){

    if(is.na(Z[i,j])==T){

      Z[i,j] <- rnorm(1)

    }
  }
}

Y <- ifelse(is.na(bin2)==T, 0, 1)*1000 + matrix(rnorm(n*p, 0 , sd=sqrt(tau2)), ncol=p, nrow=n)

sum(is.na(bin2))/(n*p)*100 # 3.371429% missing value

# Gibbs -------------------------------------------------------------------

d = 3

set.seed(125)
U <- rmnorm(n, rep(0, d), diag(1,d))
V <- t(rmnorm(p, rep(0, d), diag(1,d)))
R <- rmnorm(n, rep(0, d), diag(1,d))
W <- t(rmnorm(p, rep(0, d), diag(1,d)))

Z <- (Z-mean(Z))/sd(Z)
V <- (V-mean(V))/sd(V)
U <- (U-mean(U))/sd(U)
Y <- (Y-mean(Y))/sd(Y)
W <- (W-mean(W))/sd(W)
R <- (R-mean(R))/sd(R)
mU <- mR <- mV <- mW <- rep(0, d)
Mh <- 1e-5
Mf <- 1e+5
set.seed(125)
nIter <- 5000
gb_sen <- gibbs(nIter, U, R, V, W, Z, Y, sigma2, tau2, "binary")


# Plots -------------------------------------------------------------------

load("gb_sen.RData")
nIter <- 5000
n=100; p=35
(row_c <- similarity_funct(n, nIter, burn = 0.5*nIter, gb_sen$partition_list, dim = 1))
(col_c <- similarity_funct(p, nIter, burn = 0.5*nIter, gb_sen$partition_list, dim = 2))

assign <- row_c[[1]][3,]
nomi <- as.data.frame(cbind(data2$name, assign))
colnames(nomi) <- c("name", "cluster")

elenco <- read.csv("elenco.csv",sep=";")[-1]

results <- merge(nomi, elenco, id="name")

results$assign <- as.factor(results$cluster)
results$Party <- as.factor(results$Party)

res <- results %>%
  group_by(Party, cluster) %>%
  summarise(n=n()) %>%
  arrange(desc(Party))

res$Party <- as.factor(res$Party)

res$Party <- factor(res$Party , levels=c("Republican", "Independent", "Democrat"))

lab <- c(1:5, "Republican", "Independent", "Democrat")

# Democrats red, Republican blue
ggplot(res, aes(y = n, axis1 = Party, axis2 = cluster)) +
  geom_alluvium(aes(fill = sort(Party)), width = 1/10) +
  geom_stratum(width = 1/10, fill = "lightgrey", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size=5) +
  scale_x_discrete(limits = c("Party", "Clusters"), expand = c(.05, .05))+
  scale_fill_manual(values = alpha(c("#F8766D", "#00BA38", "#619CFF")), name = "Party")+
  ggtitle("Cluster allocation - Senators")+
  theme_minimal()+
  theme(legend.key.size = unit(2, 'cm'),
        axis.text.y=element_blank(),
        plot.title = element_text(size=22, hjust = 0.5),
        axis.title = element_text(size=20),
        axis.text.x = element_text(size = 20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18)
  )

votes <- c(paste0("Vote0", 1:9), paste0("Vote", 10:p))
type <- rep("Motion", 35)
type[c(2,4,19,21,23,27,30,31,33,34)] <- "Nomination"
type[5] <- "Resolution"

col_res <- data.frame(cluster = as.factor(col_c[[1]][3,]), votes = votes, type = type)

ggplot(col_res, aes(y = p, axis1 = type, axis2 = cluster)) +
  geom_alluvium(aes(fill = type), width = 1/20) +
  geom_stratum(width = 1/10, fill = "lightgrey", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size=5) +
  scale_x_discrete(limits = c("Votes", "Clusters"), expand = c(.05, .05))+
  scale_color_viridis(discrete = TRUE)+
  scale_fill_viridis(discrete = TRUE, option = "D", name = "Topic") +
  ggtitle("Cluster allocation - Senators' votes")+
  theme_minimal()+
  theme(legend.key.size = unit(2, 'cm'),
        axis.text.y=element_blank(),
        plot.title = element_text(size=22, hjust = 0.5),
        axis.title = element_text(size=20),
        axis.text.x = element_text(size = 20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18)
  )

# Competitor ---------------------------------------------------------------

library(biclustermd)
set.seed(1234)
bcm <- biclustermd(data = bin2, col_clusters = 2, row_clusters = 3,
                   miss_val = mean(bin2, na.rm = TRUE), miss_val_sd = 1,
                   col_min_num = 1, row_min_num = 10,
                   col_num_to_move = 1, row_num_to_move = 2,
                   col_shuffles = 1, row_shuffles = 2,
                   max.iter = 1000)
bcm
autoplot(bcm$Similarities, ncol = 3) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = 0:9)

autoplot(bcm$SSE) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = 0:9)

autoplot(bcm, reorder = TRUE, transform_colors = TRUE, c = 1/15) +
  scale_fill_viridis_c(na.value = 'white') +
  labs(x = "Senators",
       y = "Votes",
       fill = "Average Votes")

g_col <- bcm$P[,1]; g_col[g_col==0] <- 2

arandi(g_col, col_res$cluster)
table(g_col, col_res$cluster)

g_row <- bcm$Q[,1]; g_row[bcm$Q[,2]==1] <- 3; g_row[bcm$Q[,3]==1] <- 1; g_row[bcm$Q[,1]==1] <- 2

table(g_row)
table(assign)
arandi(g_row, assign)
table(g_row, assign)
nomi2 <- as.data.frame(cbind(data2$name, g_row))
colnames(nomi2) <- c("name", "cluster")

results2 <- merge(nomi2, elenco, id="name")

results2$assign <- as.factor(results2$cluster)
results2$Party <- as.factor(results2$Party)

res2 <- results2 %>%
  group_by(Party, cluster) %>%
  summarise(n=n()) %>%
  arrange(desc(Party))

res2$Party <- as.factor(res2$Party)

res2$Party <- factor(res2$Party , levels=c("Republican", "Independent", "Democrat"))

lab <- c(1:5, "Republican", "Independent", "Democrat")

ggplot(res2, aes(y = n, axis1 = Party, axis2 = cluster)) +
  geom_alluvium(aes(fill = sort(Party)), width = 1/10) +
  geom_stratum(width = 1/10, fill = "lightgrey", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size=5) +
  scale_x_discrete(limits = c("Party", "Clusters"), expand = c(.05, .05))+
  scale_fill_manual(values = alpha(c("#F8766D", "#00BA38", "#619CFF")), name = "Party")+
  ggtitle("Cluster allocation - Senators - Competitor")+
  theme_minimal()+
  theme(legend.key.size = unit(2, 'cm'),
        axis.text.y=element_blank(),
        plot.title = element_text(size=22, hjust = 0.5),
        axis.title = element_text(size=20),
        axis.text.x = element_text(size = 20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18)
  )

votes <- c(paste0("Vote0", 1:9), paste0("Vote", 10:p))
col_res2 <- data.frame(cluster = as.factor(g_col),
                       votes = votes, type=type)

ggplot(col_res2, aes(y = p, axis1 = type, axis2 = cluster)) +
  geom_alluvium(aes(fill = type), width = 1/20) +
  geom_stratum(width = 1/10, fill = "lightgrey", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size=5) +
  scale_x_discrete(limits = c("Votes", "Clusters"), expand = c(.05, .05))+
  scale_color_viridis(discrete = TRUE)+
  scale_fill_viridis(discrete = TRUE, option = "D", name = "Topic") +
  ggtitle("Cluster allocation - Senators' votes - Competitor")+
  theme_minimal()+
  theme(legend.key.size = unit(2, 'cm'),
        axis.text.y=element_blank(),
        plot.title = element_text(size=22, hjust = 0.5),
        axis.title = element_text(size=20),
        axis.text.x = element_text(size = 20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18)
  )
