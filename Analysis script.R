#### Environment setting ####

rm(list=ls())
library(car)
library(psych)
library(bootnet)
library(qgraph)
library(pcalg)
library(graph)
library(Rgraphviz)
set.seed(108)

#### Loading and processing the data ####

data  <-  read.table("Scaleultimo.csv", header<-TRUE, sep<-";", dec = ",")
data[, 2] <- factor(data[, 2], levels=c("0","1"),labels=c("Female","Male"))
data$Gender_Male <- as.numeric(data$Gender == "Male")

data$Years_of_scholarization <- data$Annidiscolarizzazione

data_backup <- data
data <- data[,c(8:40, 42:45, 3)]

colnames(data) <- c("Self blame", "Acceptance", "Rumination", "Positive\nrefocusing", "Refocus\non planning", 
                    "Positive\nreappraisal", "Putting into\nperspective", "Catastrophizing", "Other blame", "Task\noriented", 
                    "Emotion\noriented", "Treat oneself\noriented", "Contact friend\noriented", "Machiavellianism", "Psychopathy", "Narcissism", 
                    "Perception\nof self", "Perception\nof future", "Structured\nstyle", "Social\ncompetence", 
                    "Family\ncohesion", "Social\nresources", "Perceiving\nemotion", "Use of\nemotion", "Understanding\nemotion", "Managing\nemotion", 
                    "Social\nmanagement", "Covert\nnarcissism", "Extraversion", "Agreeableness", "Conscientiousness", "Emotional\nstability", "Openness", 
                    "Cognitive", "Motivational", "Gender \n(male)", "Years of \nscholarization", "Age")

#### Unregularized network ####

cortable <- cor(data)

results <- estimateNetwork(data, default="cor")
plot(results, layout="spring", groups=list(CERQ=c(1:9), CISS=c(10:13), DD=c(14:16), RSA=c(17:22), SREIS=c(23:27), NVQ=c(28), SIMP=c(29:33), SOCS=c(34:35), `Socio-biogr.`=c(36:38)), palette="pastel", legend=TRUE, legend.cex=.8, GLratio=3, 
     mar=c(1,1,1,1), vsize=7, vsize2=3.5, shape="ellipse", label.cex=1.7)

# Partial correlations
results <- estimateNetwork(data, default="pcor")
plot(results, layout="spring", groups=list(CERQ=c(1:9), CISS=c(10:13), DD=c(14:16), RSA=c(17:22), SREIS=c(23:27), NVQ=c(28), SIMP=c(29:33), SOCS=c(34:35), `Socio-biogr.`=c(36:38)), palette="pastel", legend=TRUE, legend.cex=.8, GLratio=3, 
     mar=c(1,1,1,1), vsize=7, vsize2=3.5, shape="ellipse", label.cex=1.7)

#### LASSO-regularized network ####

results <- estimateNetwork(data, default="EBICglasso", corMethod = "cor_auto")

# removal of non-significant correlations

matrix <- results$results$optwi
rownames(matrix) <- colnames(data)
colnames(matrix) <- colnames(data)
diag(matrix) <- 1
results$results$optnet[psych::corr.p(r=matrix, n=nrow(data), adj="none")$p > .05] <- 0
results$results$optwi[psych::corr.p(r=matrix, n=nrow(data), adj="none")$p > .05] <- 0
results$graph[psych::corr.p(r=matrix, n=nrow(data), adj="none")$p > .05] <- 0

plot(results, layout="spring", groups=list(CERQ=c(1:9), CISS=c(10:13), DD=c(14:16), RSA=c(17:22), SREIS=c(23:27), NVQ=c(28), SIMP=c(29:33), SOCS=c(34:35), `Socio-biogr.`=c(36:38)), palette="pastel", legend=TRUE, legend.cex=.8, GLratio=3, 
     mar=c(1,1,1,1), vsize=7, vsize2=3.5, shape="ellipse", label.cex=1.7)

centralityPlotMod(results, include=c("Strength", "Betweenness", "ExpectedInfluence"), labels=c("CERQ Self blame", "CERQ Acceptance", "CERQ Rumination", "CERQ Positive refocusing", "Refocus on planning", 
                  "CERQ Positive reappraisal", "CERQ Putting into perspective", "CERQ Catastrophizing", "CERQ Other blame", "CISS Task oriented", 
                  "CISS Emotion oriented", "CISS Treat oneself oriented", "CISS Contact friend oriented", "DD Machiavellianism", "DD Psychopathy", "DD Narcissism", 
                  "RSA Perception of self", "RSA Perception of future", "RSA Structured style", "RSA Social competence", 
                  "RSA Family cohesion", "RSA Social resources", "SREIS Perceiving emotion", "SREIS Use of emotion", "SREIS Understanding emotion", "SREIS Managing emotion", 
                  "SREIS Social management", "NVQ Covert narcissism", "SIMP Extraversion", "SIMP Agreeableness", "SIMP Conscientiousness", "SIMP Emotional stability", "SIMP Openness", 
                  "SOCS Cognitive", "SOCS Motivational", "Gender (male)", "Years of scholarization", "Age"))

#### Causal discovery ####

# Editing the pcalgo package function to make the graph more readable

newplotPC <- function (x, y, ...)
{
  .local <- function (x, y, main = NULL, zvalue.lwd = FALSE, 
                      lwd.max = 7, labels = NULL, ...) 
  {
    
    attrs <- nodeAttrs <- list()
    p <- numNodes(G <- x@graph)
    if (!is.null(labels)) {
      attrs$node <- list(shape = "ellipse", fixedsize = FALSE, height=3, cex=3, fontsize=40)
      attrs$edge <- list(arrowsize = 2, labeldistance=3)
      attrs$graph <- list(layout="dot")
      nodeAttrs <- makeNodeAttrs(G, label=c("Self blame", "Acceptance", "Rumination", "Positive \\\nrefocusing", "Refocus \\\non planning", 
                                            "Positive \\\nreappraisal", "Putting into \\\nperspective", "Catastrophizing", "Other blame", "Task \\\noriented", 
                                            "Emotion \\\noriented", "Treat oneself \\\noriented", "Contact friend \\\noriented", "Machiavellianism", "Psychopathy", "Narcissism", 
                                            "Perception \\\nof self", "Perception \\\nof future", "Structured \\\nstyle", "Social \\\ncompetence", 
                                            "Family \\\ncohesion", "Social \\\nresources", "Perceiving \\\nemotion", "Use of \\\nemotion", "Understanding \\\nemotion", "Managing \\\nemotion", 
                                            "Social \\\nmanagement", "Covert \\\nnarcissism", "Extraversion", "Agreeableness", "Conscientiousness", "Emotional \\\nstability", "Openness", 
                                            "Cognitive", "Motivational", "Gender \\\n(male)", "Years of \\\nscholarization", "Age"), shape="ellipse", fillcolor=c(rep("#FFC5D0", 9),
                                rep("#F2CDB2", 4),
                                rep("#D4D8A7", 3),
                                rep("#AFE0B9", 6),
                                rep("#99E2D8", 5),
                                rep("#ABDBF3", 1),
                                rep("#D5D0FC", 5),
                                rep("#F5C6EE", 2),
                                rep("#F5BBBB", 3)))
    }
    if (zvalue.lwd && numEdges(G) != 0) {
      lwd.mat <- if (is.matrix(Z <- x@zMin) && all(dim(Z) == 
                                                   p)) 
        Z
      else qnorm(x@pMax/2, lower.tail = FALSE)
      lwd.mat <- lwd.max * lwd.mat/max(lwd.mat)
      z <- Rgraphviz::agopen(G, name = "lwdGraph", 
                             nodeAttrs = nodeAttrs, attrs = attrs)
      for (i in seq_along(z@AgEdge)) {
        z@AgEdge[[i]]@lwd <- lwd.mat[as.integer(z@AgEdge[[i]]@head), 
                                     as.integer(z@AgEdge[[i]]@tail)]
        z@AgEdge[[i]]@col <- as.integer(z@AgEdge[[i]]@head) < 0 + 2
      }
      Rgraphviz::plot(z, main = main, ...)
    }
    else {
      Rgraphviz::plot(G, nodeAttrs = nodeAttrs, main = main, 
                      attrs = attrs, ...)
    }
  }
  .local(x, y, ...)
}


# Causal discovery

suffStat <- list(C = cor(data, use="pair"), n = nrow(data))


pc.gmG <- pc(suffStat, indepTest = gaussCItest,
             p = ncol(data), alpha = .01, conservative=FALSE, maj=FALSE, solve= TRUE)

newplotPC(pc.gmG,  labels=c("Self blame", "Acceptance", "Rumination", "Positive refocusing", "Refocus on planning", 
                            "Positive reappraisal", "Putting into perspective", "Catastrophizing", "Other blame", "Task oriented", 
                            "Emotion oriented", "Treat oneself oriented", "Contact friend oriented", "Machiavellianism", "Psychopathy", "Narcissism", 
                            "Perception of self", "Perception of future", "Structured style", "Social competence", 
                            "Family cohesion", "Social resources", "Perceiving emotion", "Use of emotion", "Understanding emotion", "Managing emotion", 
                            "Social management", "Covert narcissism", "Extraversion", "Agreeableness", "Conscientiousness", "Emotional stability", "Openness", 
                            "Cognitive", "Motivational", "Gender (male)", "Years of scholarization", "Age"))


# Bootstrap for edge weight and centrality measures
set.seed(108)
boot <- bootnet(estimateNetwork(data, default="EBICglasso", corMethod = "cor_auto"), nBoots=3000, nCores=8)
plot(boot, labels = FALSE, order = "sample")

boot2 <- bootnet(estimateNetwork(data, default="EBICglasso", corMethod = "cor_auto"), nBoots=3000, type = "case", statistics = c("strength", "betweenness", "expectedInfluence"), nCores=8)
plot(boot2, statistics = c("strength", "betweenness", "expectedInfluence"))

corStability(boot2)
