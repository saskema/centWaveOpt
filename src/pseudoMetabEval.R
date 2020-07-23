## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see https://www.gnu.org/licenses/.


## Load necessary packages
library("tidyverse")

## Define function for univariate statistics
uniStats <- function(features, names, classes, class.1, class.2) {
    
    features.class.1 <- as.matrix(features[, classes == class.1])
    features.class.2 <- as.matrix(features[, classes == class.2])
    
    summary <- as.data.frame(matrix(data = names, nrow = length(features[, 1]),
                                    dimnames = list(NULL, c("feature"))))
    
    for (i in 1:length(features[, 1])) {
        
        summary$foldchange[i] <- mean(features.class.2[i, ])/mean(features.class.1[i,])
        
    }
    for (i in 1:length(features[, 1])) {
        
        t <- t.test(features.class.1[i, ], features.class.2[i, ])
        summary$tstat[i] <- t$statistic
        summary$pvalue[i] <- t$p.value
        
    }
    
    return(summary)
    
}

## Import table
peaklist <- read.csv("features.csv", header = TRUE, dec = ".", sep = ",")
peaklist <- column_to_rownames(peaklist, "X")

Classes <- c("Blank", "Blank", "Blank", "Blank", "Blank", "Mix High", "Mix High", 
    "Mix High", "Mix High", "Mix High", "Mix Low", "Mix Low", "Mix Low", "Mix Low", 
    "Mix Low")
class.levels <- c("Blank", "Mix High", "Mix Low")

## Extract feature abundances and replace 0 by lowest measured abundance (acc. to
## Wehrens 2016)
matrix.unfiltered <- peaklist[, (10 + length(class.levels)):(length(peaklist) - 3)]

## Filtering by significance in volcano plots
uniStats.1.2 <- uniStats(matrix.unfiltered, rownames(peaklist), Classes, class.levels[1],
                         class.levels[2])
uniStats.1.3 <- uniStats(matrix.unfiltered, rownames(peaklist), Classes, class.levels[1],
                         class.levels[3])
uniStats.2.3 <- uniStats(matrix.unfiltered, rownames(peaklist), Classes, class.levels[2],
                         class.levels[3])
uniStats.1.2$Significant <- ifelse(uniStats.1.2$pvalue < 0.001 &
                                   (uniStats.1.2$foldchange > 2 |
                                    uniStats.1.2$foldchange == Inf),
                                   "TRUE", "FALSE")
uniStats.1.3$Significant <- ifelse(uniStats.1.3$pvalue < 0.001 &
                                   (uniStats.1.3$foldchange > 2 |
                                    uniStats.1.3$foldchange == Inf),
                                   "TRUE", "FALSE")
uniStats.2.3$Significant <- ifelse(uniStats.2.3$pvalue < 0.001 &
                                   (uniStats.2.3$foldchange > 2 |
                                    uniStats.2.3$foldchange == Inf),
                                   "TRUE", "FALSE")

uniStats.significant <- list(uniStats_1 = uniStats.1.2$feature[uniStats.1.2$Significant == TRUE],
                             uniStats_2 = uniStats.1.3$feature[uniStats.1.3$Significant == TRUE],
                             uniStats_3 = uniStats.2.3$feature[uniStats.2.3$Significant == TRUE])
uniStats.significant$all <- purrr::reduce(uniStats.significant, union)


lod <- min(matrix.unfiltered[matrix.unfiltered != 0])
matrix.unfiltered[matrix.unfiltered == 0] <- lod
matrix.unfiltered <- log10(matrix.unfiltered)

matrix.filtered <- column_to_rownames(filter(rownames_to_column(as.data.frame(matrix.unfiltered)),
                                             peaklist$Low == 5, peaklist$rt > 60, peaklist$rt < 720))
matrix.filtered <- t(matrix.filtered)


## Perform PCA
pca <- prcomp(matrix.filtered, scale. = FALSE, center = TRUE)
sum <- summary(pca)
peaklist.filtered <- column_to_rownames(filter(rownames_to_column(as.data.frame(peaklist)),
                                               peaklist$Low == 5, peaklist$rt > 60, peaklist$rt < 720))
loadings <- as.data.frame(pca$rotation[, 1:2])
loadings$Category <- peaklist.filtered$category
loadings$Significant <- ifelse(colnames(matrix.filtered) %in% uniStats.significant$all, "Yes", "No")


## Plotting Results
colours_classes <- c(Blank = "#619CFF", `Mix High` = "#F8766D", `Mix Low` = "#00BA38")
colours_features <- c(`Parent Compound` = "#F8766D", Isotope = "#00BA38", Adduct = "#DAA524",
                      Artifact = "#C77CFF", Unknown = "#619CFF")

## Plot PCA Scores
p.scores <- ggplot(data = as.data.frame(pca$x[, 1:2]), aes(x = PC1, y = PC2, col = Classes)) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) + 
    scale_color_manual(values = colours_classes) +
    xlab(paste0("Principal Component 1", 
                " ", "(", round(sum$importance[2, 1] * 100, 0), "%)")) +
    ylab(paste0("Principal Component 2", 
                " ", "(", round(sum$importance[2, 2] * 100, 0), "%)")) +
    theme_bw()

## Plot PCA Loadings
p.loadings <- ggplot(data = loadings, aes(x = PC1, y = PC2, col = Category, shape = Significant)) + 
    geom_point(size = 3) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) + 
    scale_color_manual(values = colours_features) +
    xlab(paste0("Principal Component 1", 
                " ", "(", round(sum$importance[2, 1] * 100, 0), "%)")) +
    ylab(paste0("Principal Component 2", 
                " ", "(", round(sum$importance[2, 2] * 100, 0), "%)")) +
    theme_bw()

ggsave("Scores.pdf", p.scores, width = 6, height = 4, device = "pdf")
ggsave("Loadings.pdf", p.loadings, width = 6, height = 4, device = "pdf")
