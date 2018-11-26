library(ggplot2)
library(tidyr)

#Enterococcus
tab_data_e <- read.csv("~/Desktop/FYP/resultsm/enterococcus/ANI/ANIm_percentage_identity.tab", sep = "\t", stringsAsFactors = F, )
row.names(tab_data_e)<-tab_data_e$X
tab_data_e$X <- NULL
str(tab_data_e)
values_e <- c()
indexes_of_interest_e <- c() 
for (col in 1:ncol(tab_data_e)){
  name <- colnames(tab_data_e[col])
  if (startsWith(name, "X.NODE")){
    #print(name)
    if (is.na(tab_data_e[col-1, col])){
      next
    }
    if (is.na(tab_data_e[col+1, col])){
      next
    }
    indexes_of_interest_e <- c(indexes_of_interest_e,
                             col, col-1, col+1)
    values_e <- c(values_e,tab_data_e[col-1, col], tab_data_e[col+1, col] )
    if (tab_data_e[col-1, col] < .99){
       print("hit above")
     } else if (tab_data_e[col + 1, col] < .99){
      print("hit below")
     }
  }
    
}

small_matrix_e <- tab_data_e[unique(indexes_of_interest_e),
                         unique(indexes_of_interest_e)]

hist(values_e)
tall <- tab_data_e %>% gather(key, value)
ggplot(tall, aes(x=X, y=key)) + geom_tile()

#install.packages("gplots")
#library(gplots)
#heatmap.2(as.matrix(tab_data_e),
#          breaks=c(.90, .99, .995, 1))

#heatmap.2(as.matrix(small_matrix_e),
#          breaks=c(.90, .99, .995, 1),
#          Rowv = row.names(small_matrix_e), Colv = colnames(small_matrix_e))
        
heatmap.2(as.matrix(small_matrix_e),
          density.info = "none",
          trace  = "column",
            tracecol = "#555555",
          key.title = "Percentage similarity",
            
          breaks=c(.98,.985,.99,.993,.999, 1), 
          col = c("#111111","white","#333333","#555555","#999999"),
          Rowv = row.names(small_matrix_e), Colv = colnames(small_matrix_e))

heat.colors

#Morganii 

tab_data_m <- read.csv("~/Desktop/FYP/resultsm/morganii/ANI/ANIm_percentage_identity.tab", sep = "\t", stringsAsFactors = F, )
row.names(tab_data_m)<-tab_data_m$X
tab_data_m$X <- NULL
str(tab_data_m)
values_m <- c()
indexes_of_interest_m <- c() 
for (col in 1:ncol(tab_data_m)){
  name <- colnames(tab_data_m[col])
  if (startsWith(name, "X.NODE")){
    #print(name)
    ABOVE_HIT = FALSE
    BELOW_HIT = F
    values_m <- c(values_m,tab_data_m[col-1, col], tab_data_m[col+1, col] )
    if (col != 1){
      indexes_of_interest_m <- c(indexes_of_interest_m,
                                 col, col-1)
      if (tab_data_m[col-1, col] < .99){
        ABOVE_HIT = T
      }
    }
    if (col != ncol(tab_data_m)){
      indexes_of_interest_m <- c(indexes_of_interest_m,col, col+1) 
      if (tab_data_m[col + 1, col] < .99){
        BELOW_HIT = T
      }
    }
    if (BELOW_HIT & ABOVE_HIT){
      print("novel sequence")
      print(colnames(tab_data_m[col]))
    }
  }
}

small_matrix_m <- tab_data_m[unique(indexes_of_interest_m),
                         unique(indexes_of_interest_m)]
# reduce to single triangle
get_upper_tri <- function(cormat){
  cormat[upper.tri(cormat)]<- NA
  return(cormat)
}
small_matrix_m <- get_upper_tri(small_matrix_m)

#hist(values_m)
tall <- small_matrix_m %>% mutate(X=colnames(small_matrix_m)) %>% gather(key, value, -X)
#reorder x
tall$X <- factor(tall$X, levels=colnames(small_matrix_m))
# reorder y
tall$key <- factor(tall$key, levels=colnames(small_matrix_m))
#plot
tall$riboseed <- ifelse(tall$key == tall$X &
                          tall$X %in% c("X.NODE_2_length_813131_cov_12.3972.RC_",
                                     "X.NODE_4_length_288872_cov_12.2011.RC_",
                                     "X.NODE_13_length_84248_cov_13.2788.RC_"),
                      T, F) 

 ggplot(tall, aes(x=X, y=key, fill=value, color=riboseed)) +
  geom_tile(size=.8) +
  theme(
    axis.text.x = element_text(angle=45, hjust = 1),
    legend.position=c(.1, .6)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = .97, limit = c(.95,1), space = "Lab", 
                       name="Percentage\nSimilarity", na.value = "white") +
   scale_color_manual(values=c("white", "black"), guide="none")+
   labs(x="Sequence",
       y="Sequence") +
   scale_y_discrete(position = "right") 

     
  

#install.packages("gplots")
#library(gplots)
#heatmap.2(as.matrix(tab_data_m),
#          breaks=c(.90, .99, .995, 1))
#
#heatmap.2(as.matrix(small_matrix_m),
#          breaks=c(.90, .99, .995, 1),
#          Rowv = row.names(small_matrix_m), Colv = colnames(small_matrix_m))

# heatmap.2(as.matrix(small_matrix_m),
#           density.info = "none",
#           trace  = "none",
#           dendrogram = "none",
#             tracecol = "#555555",
#           key.title = "Percentage similarity",
#           breaks=c(0, .97,.99, 1),
#           col = c("white","#FC0D1B","#0b24fa"),
#           Rowv = row.names(small_matrix_m), Colv = colnames(small_matrix_m))
heat.colors
