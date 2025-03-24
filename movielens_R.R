rm(list=ls())

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(dslabs)
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
library(ggalluvial)

source("gibbs_movie.R")

# Data --------------------------------------------------------------------

data(movielens)
head(movielens)

data_full <- movielens %>%
  select(userId, movieId, rating) %>%
  pivot_wider(values_from = "rating", names_from="movieId")

data_full <- data_full[,-1] # remove userId

count <- c()
for(i in 1:nrow(data_full)){

  count[i] <- sum(!is.na(data_full[i,]))

}

df1 <- data_full[count>350,]

count2 <- c()
for(j in 1:ncol(df1)){

  count2[j] <- sum(!is.na(df1[,j]))

}

df2 <- df1[,count2>45]
sum(is.na(df2))/(nrow(df2)*ncol(df2))*100 # 18%
plot(as.matrix(df2), breaks= c(0,cutoff,5.25),
     main="Movielens data", xlab="Users", ylab="Movies")

dim(df2) # 60x28

index_title <- colnames(df2)

title <- c()
genre <- c()
for(el in index_title){

  title <- c(title, movielens$title[which(movielens$movieId==el)][1])
  genre <- c(genre, as.character(movielens$genres[which(movielens$movieId==el)][1]))

}

title[1] <- "Seven"
title[6] <- "The Fugitive"
title[8] <- "The Silence of the Lambs"
title[11] <- "Men in Black"
title[9] <- "The Shawshank Redemption"
title[10] <- "Star Wars: Ep. VI"
title[11] <- "Men in Black"
title[12] <- "The Sixth Sense"
title[14] <- "Star Wars: Ep. IV"
title[15] <- "The Godfather"
title[19] <- "Star Wars: Ep. V"
title[20] <- "Raiders of the Lost Ark"
title[22] <- "The Terminator"
title[26] <- "The Matrix"

library(colorBlindness)
par(pty="s")
mybar <- barplot(table(unlist(as.vector(df2)[!is.na(as.vector(df2))])),
          ylim=c(0,470), col="grey",
          xlab="Ratings", main="Counts per ratings", ylab="Counts"
        )
text(mybar, table(unlist(as.vector(df2)[!is.na(as.vector(df2))]))+22 ,
                   table(unlist(as.vector(df2)[!is.na(as.vector(df2))])) ,cex=1)

library(ggplot2)
library(dplyr)
library(tidyr)

ratings <- unlist(as.vector(df2))
ratings <- ratings[!is.na(ratings)]

ratings_count <- as.data.frame(table(ratings))

ggplot(ratings_count, aes(x = ratings, y = Freq)) +
  geom_bar(stat = "identity", fill = "grey", color = "black") +
  geom_text(aes(label = Freq), vjust = -0.5, size = 5) +
  ylim(0, 470) +
  labs(x = "Ratings", y = "Counts") +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
      legend.key.size = unit(2, 'cm'),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(size=22, hjust = 0.5),
      axis.title = element_text(size=20),
      axis.text = element_text(size = 20),
      legend.title=element_text(size=20),
      legend.text=element_text(size=18))


df2 <- as.matrix(df2)
# cutoff
cutoff <- seq(0.75, 5, 0.5)

n <- nrow(df2)
p <- ncol(df2)
rownames(df2) <- c(1:n)
colnames(df2) <- c(1:p)
par(mar=c(5.1,4.1,4.1,4.1))
plot(df2, xaxt="n", yaxt="n", main="Movielens ratings",
     xlab="Movies", ylab="Users", breaks=c(0.25,cutoff, 5.25))

set.seed(125)
sigma2 <- 0.1
tau2 <- 1.5
Z <- df2*1000 + matrix(rnorm(n*p, 0 , sd=sqrt(sigma2)), ncol=p, nrow=n)

for(i in 1:n){
  for(j in 1:p){

    if(is.na(Z[i,j])==T){

      Z[i,j] <- rnorm(1)

    }
  }
}

Y <- ifelse(is.na(df2)==T, 0, 1)*1000 + matrix(rnorm(n*p, 0 , sd=sqrt(tau2)), ncol=p, nrow=n)
sum(is.na(df2))/(n*p)*100 # 17.97% missing value

Z <- (Z-mean(Z))/sd(Z)
Y <- (Y-mean(Y))/sd(Y)
plot(Z, breaks=c(-2, scale(seq(0.5, 5, by=0.5)), 2))

sc_cutoff <- scale(seq(0.5, 5, by=0.5))+0.1651446

# Lpml --------------------------------------------------------------------
#
# D <- 2:20
#
# gb_lpml <- list()
# for(d in D){
#
#   set.seed(125)
#   U <- rmnorm(n, rep(0, d), diag(1,d))
#   V <- t(rmnorm(p, rep(0, d), diag(1,d)))
#   R <- rmnorm(n, rep(0, d), diag(1,d))
#   W <- t(rmnorm(p, rep(0, d), diag(1,d)))
#
#   V <- (V-mean(V))/sd(V)
#   U <- (U-mean(U))/sd(U)
#   W <- (W-mean(W))/sd(W)
#   R <- (R-mean(R))/sd(R)
#
#   Mh <- 1
#   Mf <- 1
#   nIter <- 1000
#   set.seed(125)
#
#
#   mU <- mR <- mV <- mW <- rep(0, d)
#
#   gb_lpml[d] <- list(gibbs_cpo(nIter, U, R, V, W, Z, Y, sigma2, tau2, "ordinal"))
#
# }
#
# lpml <- c()
# for(d in 2:20){
#   lpml[d-1] <- gb_lpml[[d-1]]$lpml
# }
#
# plot(lpml[complete.cases(lpml)], type="l", xaxt="n", main="Latent dimension d", ylab = "LPML")
# axis(1, at=1:19, labels=D)
# abline(v=which(lpml[complete.cases(lpml)]==max(lpml[complete.cases(lpml)])), col=2)
load("lpml_movie.RData")
par(pty="s", mar=c(5,6,4,2)+0.1)
plot(lpml, xaxt="n", xlab="d",
     ylab = "LPML", pch=20, main="", cex=2,
     cex.lab=2.5, cex.axis=1.5, cex.sub=1.5, cex.main=2)
abline(v=1:19, col="grey", lwd=1, cex=2)
points(lpml, pch=16, cex=2)
axis(1, at=1:19, labels=D, las=1, cex.axis=1.5)

# Gibbs -------------------------------------------------------------------

d = 2

set.seed(125)
U <- rmnorm(n, rep(0, d), diag(1,d))
V <- t(rmnorm(p, rep(0, d), diag(1,d)))
R <- rmnorm(n, rep(0, d), diag(1,d))
W <- t(rmnorm(p, rep(0, d), diag(1,d)))

V <- (V-mean(V))/sd(V)
U <- (U-mean(U))/sd(U)
W <- (W-mean(W))/sd(W)
R <- (R-mean(R))/sd(R)
mU <- mR <- mV <- mW <- rep(0, d)
Mh <- 1e-5
Mf <- 1e+50
set.seed(125)
nIter <- 5000
data <- df2
sigm2 <- 0.1
tau2 <- 1.5
gb_movie <- gibbs(nIter, U, R, V, W, Z, Y, sigma2, tau2, "ordinal")

# Plots -------------------------------------------------------------------

(row_c <- similarity_funct(n, nIter, burn = 0.5*nIter, gb_movie$partition_list, dim = 1))
(col_c <- similarity_funct(p, nIter, burn = 0.5*nIter, gb_movie$partition_list, dim = 2))

assign <- row_c[[1]][3,]

res <- data.frame(Id = 1:n, cluster=assign)

ggplot(res, aes(y = Id, axis1 = 1:n, axis2 = as.factor(cluster))) +
  geom_alluvium(aes(fill = as.factor(cluster)), width = 1/20) +
  geom_stratum(width = 1/10, fill = "lightgrey", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Users", "Clusters"), expand = c(.005, .05))+
  scale_color_viridis(discrete = TRUE,  option = "C", name = "Cluster")+
  scale_fill_viridis(discrete = TRUE, option = "C", name = "Cluster") +
  ggtitle("Cluster allocation - Users")+
  ggrepel::geom_text_repel(
    aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), NA)),
    stat = "stratum", size = 2, direction = "y", nudge_x = -.12)+
  theme(legend.key.size = unit(2, 'cm'), axis.text.y=element_blank(),
        plot.title = element_text(size=22, hjust = 0.5),
        axis.title = element_text(size=20),
        axis.text.x = element_text(size = 20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18))

cluster <- as.factor(col_c[[1]][3,])
col_res <- data.frame(cluster = cluster, title = title)

df_genre <- cbind(title,genre, cluster)

ggplot(col_res, aes(y = p, axis1 = title, axis2 = cluster)) +
  geom_alluvium(aes(fill = cluster), width = 1/20) +
  geom_stratum(width = 1/10, fill = "lightgrey", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Movies", "Clusters"), expand = c(.05, .05))+
  scale_color_viridis(discrete = TRUE)+
  scale_fill_viridis(discrete = TRUE, option = "D", name = "Genre") +
  ggtitle("Cluster allocation - Movies")+
  theme(legend.key.size = unit(2, 'cm'), axis.text.y=element_blank(),
        plot.title = element_text(size=22, hjust = 0.5), axis.title = element_text(size=20),
        axis.text.x = element_text(size = 20), legend.title=element_text(size=20),
        legend.text=element_text(size=18))


# Competitor ---------------------------------------------------------------

library(biclustermd)
set.seed(1234)
bcm <- biclustermd(data = data, col_clusters = 4, row_clusters = 6,
                   miss_val = mean(data, na.rm = TRUE), miss_val_sd = 1,
                   col_min_num = 2, row_min_num = 1,
                   col_num_to_move = 2, row_num_to_move = 1,
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
  labs(x = "Movies",
       y = "Users",
       fill = "Average Ratings")

g_row <- rowSums(bcm$Q%*%diag(c(1,4,3,2,6,5)))

table(g_row); table(assign)
table(g_row, assign)

res <- data.frame(Id = 1:n, cluster=g_row)

ggplot(res, aes(y = Id, axis1 = 1:n, axis2 = as.factor(cluster))) +
  geom_alluvium(aes(fill = as.factor(cluster)), width = 1/20) +
  geom_stratum(width = 1/10, fill = "lightgrey", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Users", "Clusters"), expand = c(.005, .05))+
  scale_color_viridis(discrete = TRUE,  option = "C", name = "Cluster")+
  scale_fill_viridis(discrete = TRUE, option = "C", name = "Cluster") +
  ggtitle("Cluster allocation - Users")+
  ggrepel::geom_text_repel(
    aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), NA)),
    stat = "stratum", size = 2, direction = "y", nudge_x = -.12)+
  theme(legend.key.size = unit(2, 'cm'), axis.text.y=element_blank(),
        plot.title = element_text(size=22, hjust = 0.5),
        axis.title = element_text(size=20),
        axis.text.x = element_text(size = 20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18))
#ggsave("row_movie_competitor.pdf", width = 18, height = 12, scale =0.8)

g_col <- rowSums(bcm$P%*%diag(c(4,2,1,3)))

table(g_col); table(cluster)
table(g_col, cluster)

col_res <- data.frame(cluster = as.factor(g_col), title = as.factor(title))

ggplot(col_res, aes(y = p, axis1 = title, axis2 = cluster)) +
  geom_alluvium(aes(fill = cluster), width = 1/20) +
  geom_stratum(width = 1/10, fill = "lightgrey", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Movies", "Clusters"), expand = c(.05, .05))+
  scale_color_viridis(discrete = TRUE)+
  scale_fill_viridis(discrete = TRUE, option = "D", name = "Genre") +
  ggtitle("Cluster allocation - Movies")+
  theme(legend.key.size = unit(2, 'cm'), axis.text.y=element_blank(),
        plot.title = element_text(size=22, hjust = 0.5), axis.title = element_text(size=20),
        axis.text.x = element_text(size = 20), legend.title=element_text(size=20),
        legend.text=element_text(size=18))

df_comp_movie <- cbind(title, genre, g_col)

# Radar chart -------------------------------------------------------------

# Library
library(fmsb)

dff <- matrix(NA, ncol=length(unique(col_res$cluster)), nrow=length(unique(assign)))
missperc <- matrix(NA, ncol=length(unique(col_res$cluster)), nrow=length(unique(assign)))
for(clus in unique(col_res$cluster)){

  for(cluss in unique(assign)){

    index <- which(assign==cluss)

    if(length(index)>1){

      dff[cluss,as.numeric(clus)] <- mean(rowMeans(df2[index,which(col_res$cluster==clus)], na.rm = T))

      missperc[cluss,as.numeric(clus)] <- sum(is.na(df2[index,which(col_res$cluster==clus)]))

    }else{
      dff[cluss,as.numeric(clus)] <- mean(df2[index,which(col_res$cluster==clus)], na.rm = T)
      missperc[cluss,as.numeric(clus)] <- sum(is.na(df2[index,which(col_res$cluster==clus)]))

    }

  }

}

dff[is.na(dff)] <- 0

dff <- as.data.frame(dff)

maxmin <- data.frame("Drama/Thriller"=c(5,0), "Adventure/Action"=c(5,0),
                     "Journey"=c(5,0), "Satire"=c(5,0))
colnames(maxmin) <- c("Drama/Thriller", "Adventure/Action",
                      "Journey", "Satire")
colnames(dff) <- c("Drama/Thriller", "Adventure/Action",
                   "Journey", "Satire")

dff <- rbind(maxmin, dff)

rownames(dff) <- c("Max", "Min", "User.1", "User.2",
                   "User.3", "User.4", "User.5", "User.6")
dff <- dff[,c("Drama/Thriller","Journey",
              "Adventure/Action","Satire")]
# Define colors and titles
colors <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF", "#E16462FF", "#FCA636FF", "#F0F921FF")#c("#00AFBB", "#E7B800", "#FC4E07", "")
titles <- c("User.1", "User.2", "User.3", "User.4", "User.5", "User.6")

# Reduce plot margin using par()
# Split the screen in 3 parts
op <- par(mar = c(1, 1, 1, 1))
par(mfrow = c(3,2))

create_beautiful_radarchart <- function(data, color = "#00AFBB",
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 0,
    #seg=10,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 2,
    # Customize the axis
    axislabcol = "grey",
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

# Create the radar chart
for(i in 1:6){
  create_beautiful_radarchart(
    data = dff[c(1, 2, i+2), ], caxislabels = c(0,cutoff-0.25, 5),
    color = colors[i], title = titles[i]
  )
}
par(op)

# missing

miss2 <- as.data.frame(missperc)

maxmin2 <- data.frame("Drama/Thriller"=c(50,0), "Adventure/Action"=c(55,0),
                     "Journey"=c(20,0), "Satire"=c(15,0))
colnames(maxmin2) <- c("Drama/Thriller", "Adventure/Action",
                      "Journey", "Satire")
colnames(miss2) <- c("Drama/Thriller", "Adventure/Action",
                   "Journey", "Satire")

miss2 <- rbind(maxmin2, miss2)

rownames(miss2) <- c("Max", "Min", "User.1", "User.2",
                   "User.3", "User.4", "User.5", "User.6")
miss2 <- miss2[,c("Drama/Thriller","Journey",
                  "Adventure/Action","Satire")]

# Define colors and titles
colors <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF", "#E16462FF", "#FCA636FF", "#F0F921FF")#c("#00AFBB", "#E7B800", "#FC4E07", "")
titles <- c("User.1", "User.2", "User.3", "User.4", "User.5", "User.6")

# Reduce plot margin using par()
# Split the screen in 3 parts
op <- par(mar = c(1, 1, 1, 1))
par(mfrow = c(3,2))

create_beautiful_radarchart2 <- function(data, color = "#00AFBB",
                                        vlabels = colnames(data), vlcex = 0.8,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 0,
    #seg=10,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 2,
    # Customize the axis
    axislabcol = "grey",
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

# Create the radar chart
for(i in 1:6){
  create_beautiful_radarchart2(
    data = miss2[c(1, 2, i+2), ], caxislabels = c(0,cutoff-0.25, 5),
    color = colors[i], title = titles[i]
  )
}
par(op)


