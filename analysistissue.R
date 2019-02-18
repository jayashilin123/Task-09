
library(dslabs)
library(caret)
load("tissuesGeneExpression.rda")

# plot with multidimentional scaling
d <- dist(t(e))
mds <- cmdscale(d)
plot(mds[,1], mds[,2], col=km$cluster)

plot(mds[,1], mds[,2], bg=km$cluster,pch=21,
     xlab="First dimension", ylab="Second dimension",
     cex=2) 
legend("bottomright", legend = levels(as.factor(km$cluster)),
       col=seq(along=levels(as.factor(km$cluster))), pch=15,cex=0.65)


###### Ref: Rafael Irizarry & Michael Love, High dimensional data analysis, Harvardx, edx.

library(tissuesGeneExpression)
data(tissueGeneExpression)
image(e[1:100,]) # yellow high, red low, colour represents the number, 
# each pixel/square represent each gene.

# can pick few genes to demonstrate the gene expression with help of image pixels
library(genefilter)
rv <- rowVars(e)
idx <- order(-rv)[1:40] # taking top 40 genes that vary the most, most rowvariance

heatmap(e[idx, ]) # heatmap of a subset of the 40 gene showing most variation
# so genes are clustered and samples are clustered. colours shows the measurments, 
# yellow are most expressed in that particular cluster.

library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) # green to blue, 9 color
heatmap(e[idx, ], col=hmcol) # high expression=blue, low expression green
library(gplots)
library(rafalib)
cols <- palette(brewer.pal(7, "Dark2"))[as.fumeric(tissue)]# each tissue get diff colors
cbind(colnames(e), tissue, cols)
heatmap.2(e[idx, ], labCol = tissue,
          trace= "none",
          ColSideColors = cols, # colum side colors
          col=hmcol)

# how variability makes clustering non deterministic

# Decision tree

library(dslabs)
library(caret)
load("tissuesGeneExpression.rda")
# data("tissue_gene_expression")
set.seed(1991)

fit <- with(tissue_gene_expression, 
            train(x, y, method = "rpart",
                  tuneGrid = data.frame(cp = seq(0, 0.1, 0.01))))

ggplot(fit)
confusionMatrix(fit)

plot(fit$finalModel, margin = 0.1) 
text(fit$finalModel, cex = 0.75)

# Hierarchial clustering
library(tissuesGeneExpression)
data(tissuesGeneExpression)
d <- dist(t(e))
hc <- hclust(d)
class(hc)
plot(hc, cex=0.5, label=tissue)
library(rafalib)
mypar(1,1)
myplclust(hc, cex=0.5, label=tissue,
          lab.col=as.fumeric(tissue))
abline(h=120)

cl <- cutree(hc, h=120)
table(true=tissue, cluster=cl)

# K nearest neighbour
library(dslabs)
data("tissue_gene_expression")

dim(as.matrix(tissue_gene_expression))
x <- tissue_gene_expression$x
rownames(x) <- NULL
dim(as.matrix(x))
y <- tissue_gene_expression$y
dim(as.matrix(y))

x <- as.data.frame(x)
y <- as.data.frame(y)

ks <- c(1, 3, 5, 7, 9, 11)
set.seed(1)
library(caret)
accuracy <- sapply(ks, function(k){
  test_index <- createDataPartition(y$y, times=1, p=0.5, list = FALSE)
  #sample(1:nrow(y),0.5*nrow(y))
  test_y <- y[test_index, ]
  train_y <- y[-test_index, ]
  test_x <- x[test_index, ]
  train_x <- x[-test_index, ]
  fit <- knn3(train_x, train_y, k=k)
  
  y_hat <- predict(fit, test_x, type = "class") #%>% factor(levels = levels(train_y))
  
  confusionMatrix(data = y_hat, reference = test_y)$overall["Accuracy"] 
})

accuracy
plot(ks, accuracy)

#########################
library(RColorBrewer)
library(rafalib)
d <- dist(t(e))
mds <- cmdscale(d)
plot(mds[,1], mds[,2], col=palette(brewer.pal(7, "Dark2"))[as.fumeric(tissue)])

plot(mds[,1], mds[,2], bg=palette(brewer.pal(7, "Dark2"))[as.fumeric(tissue)],pch=21,
     xlab="First dimension", ylab="Second dimension",
     cex=2) 
legend("bottomright", legend = levels(as.factor(tissue)),
       col=seq(along=levels(as.factor(tissue))), pch=15,cex=0.65)

levels(as.factor(tissue))
