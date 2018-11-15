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


## Place XCMS readable files in a folder called "QC" and enter your file
## extension!

## Load necessary packages
library("tidyverse")
library("stringr")
library("grid")
library("ggsignif")


files <- list.files(path = getwd(), recursive = TRUE, full.names = TRUE, pattern = ".csv")
filenames <- str_sub(files, start = 5, end = -5)
experiments <- c("Mix High Result A", "Mix High Result B", "Mix High Result C", "Mix High Result D",
                 "Mix Low Result A", "Mix Low Result B", "Mix Low Result C", "Mix Low Result D")

## Perform univariate statistics for features
peaklist.modified <- read.csv(files[1], header = TRUE, dec = ".", sep = ",")[, -1]
abundance.mean <- apply(peaklist.modified[, 9:18], 1, FUN = mean, na.rm = T)
abundance.sd <- apply(peaklist.modified[, 9:18], 1, FUN = sd, na.rm = T)
cv <- (peaklist.modified$abundance.sd * 100)/peaklist.modified$abundance.mean
deltart <- peaklist.modified$rtmax - peaklist.modified$rtmin
deltamz <- peaklist.modified$mzmax - peaklist.modified$mzmin
experiment <- experiments[1]
assign(filenames[1], cbind.data.frame(cv, deltart, deltamz, experiment))
summary <- get(filenames[1])

for (i in 2:length(files)) {
    
    peaklist.modified <- read.csv(files[i], header = TRUE, dec = ".", sep = ",")[,-1]
    abundance.mean <- apply(peaklist.modified[, 9:18], 1, FUN = mean, na.rm = T)
    abundance.sd <- apply(peaklist.modified[, 9:18], 1, FUN = sd, na.rm = T)
    cv <- (abundance.sd * 100)/abundance.mean
    deltart <- peaklist.modified$rtmax - peaklist.modified$rtmin
    deltamz <- peaklist.modified$mzmax - peaklist.modified$mzmin
    experiment <- experiments[i]
    assign(filenames[i], as.data.frame(cbind(cv, deltart, deltamz, experiment)))
    summary <- rbind(summary, get(filenames[i]))
}


## Summarise results
summary$cv <- as.numeric(summary$cv)
summary$deltart <- as.numeric(summary$deltart)
summary$deltart[summary$deltart == "-Inf"] <- 0
summary$deltamz <- as.numeric(summary$deltamz)
summary$deltamz[summary$deltamz == "-Inf"] <- 0
summary$experiment <- factor(summary$experiment, levels = c("Mix Low Result A", "Mix Low Result B",
                                                            "Mix Low Result C", "Mix Low Result D", "Mix High Result A",
                                                            "Mix High Result B", "Mix High Result C", "Mix High Result D"))

## Perform t-test for Mix Low experiments
t.cv.ab.low <- t.test(summary$cv[summary$experiment == "Mix Low Result A"], summary$cv[summary$experiment ==
                                                                                       "Mix Low Result B"])
t.cv.cd.low <- t.test(summary$cv[summary$experiment == "Mix Low Result C"], summary$cv[summary$experiment ==
                                                                                       "Mix Low Result D"])
t.cv.bd.low <- t.test(summary$cv[summary$experiment == "Mix Low Result B"], summary$cv[summary$experiment ==
                                                                                       "Mix Low Result D"])

t.deltart.ab.low <- t.test(summary$deltart[summary$experiment == "Mix Low Result A"],
                           summary$deltart[summary$experiment == "Mix Low Result B"])
t.deltart.cd.low <- t.test(summary$deltart[summary$experiment == "Mix Low Result C"],
                           summary$deltart[summary$experiment == "Mix Low Result D"])
t.deltart.bd.low <- t.test(summary$deltart[summary$experiment == "Mix Low Result B"],
                           summary$deltart[summary$experiment == "Mix Low Result D"])

t.deltamz.ab.low <- t.test(summary$deltamz[summary$experiment == "Mix Low Result A"],
                           summary$deltamz[summary$experiment == "Mix Low Result B"])
t.deltamz.cd.low <- t.test(summary$deltamz[summary$experiment == "Mix Low Result C"],
                           summary$deltamz[summary$experiment == "Mix Low Result D"])
t.deltamz.bd.low <- t.test(summary$deltamz[summary$experiment == "Mix Low Result B"],
                           summary$deltamz[summary$experiment == "Mix Low Result D"])

## Perform t-test for Mix High experiments
t.cv.ab.high <- t.test(summary$cv[summary$experiment == "Mix High Result A"], summary$cv[summary$experiment ==
                                                                                         "Mix High Result B"])
t.cv.cd.high <- t.test(summary$cv[summary$experiment == "Mix High Result C"], summary$cv[summary$experiment ==
                                                                                         "Mix High Result D"])
t.cv.bd.high <- t.test(summary$cv[summary$experiment == "Mix High Result B"], summary$cv[summary$experiment ==
                                                                                         "Mix High Result D"])

t.deltart.ab.high <- t.test(summary$deltart[summary$experiment == "Mix High Result A"],
                            summary$deltart[summary$experiment == "Mix High Result B"])
t.deltart.cd.high <- t.test(summary$deltart[summary$experiment == "Mix High Result C"],
                            summary$deltart[summary$experiment == "Mix High Result D"])
t.deltart.bd.high <- t.test(summary$deltart[summary$experiment == "Mix High Result B"],
                            summary$deltart[summary$experiment == "Mix High Result D"])

t.deltamz.ab.high <- t.test(summary$deltamz[summary$experiment == "Mix High Result A"],
                            summary$deltamz[summary$experiment == "Mix High Result B"])
t.deltamz.cd.high <- t.test(summary$deltamz[summary$experiment == "Mix High Result C"],
                            summary$deltamz[summary$experiment == "Mix High Result D"])
t.deltamz.bd.high <- t.test(summary$deltamz[summary$experiment == "Mix High Result B"], 
                            summary$deltamz[summary$experiment == "Mix High Result D"])

summary$cv <- log10(summary$cv)
summary$deltart <- log10(summary$deltart)
summary$deltart[summary$deltart == "-Inf"] <- 0
summary$deltamz <- log10(summary$deltamz)
summary$deltamz[summary$deltamz == "-Inf"] <- 0
summary$experiment <- factor(summary$experiment, levels = c("Mix Low Result A", "Mix Low Result B", 
                                                            "Mix Low Result C", "Mix Low Result D", 
                                                            "Mix High Result A", "Mix High Result B", 
                                                            "Mix High Result C", "Mix High Result D"))


## Plot cv of experiments in boxplot
p.cv.box <- ggplot(data = summary, aes(x = experiment, y = cv)) + geom_boxplot(aes(group = experiment)) +
    ggtitle("Coefficient of variation of peak areas in plasma samples") + xlab("Experiments") +
    ylab("Log10 of Coefficient of Variation") + 
geom_signif(comparisons = list(c("Mix Low Result A", "Mix Low Result B")), map_signif_level = TRUE, 
    test = "t.test", textsize = 6) + geom_signif(comparisons = list(c("Mix Low Result C", 
    "Mix Low Result D")), map_signif_level = TRUE, test = "t.test", textsize = 6) + 
    geom_signif(comparisons = list(c("Mix Low Result B", "Mix Low Result D")), map_signif_level = TRUE,
                test = "t.test", textsize = 6, y_position = 3) + 
geom_signif(comparisons = list(c("Mix High Result A", "Mix High Result B")), map_signif_level = TRUE, 
    test = "t.test", textsize = 6) + geom_signif(comparisons = list(c("Mix High Result C", 
    "Mix High Result D")), map_signif_level = TRUE, test = "t.test", textsize = 6) + 
    geom_signif(comparisons = list(c("Mix High Result B", "Mix High Result D")),
                map_signif_level = TRUE, test = "t.test", textsize = 6, y_position = 3) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Plot retention time range integrated by XCMS
p.deltart.box <- ggplot(data = summary, aes(x = experiment, y = deltart)) + geom_boxplot(aes(group = experiment)) + 
    ggtitle("Retention time ranges within feature groups in plasma samples") + xlab("Experiments") +
    ylab("Log10 retention time ranges") + 
geom_signif(comparisons = list(c("Mix Low Result A", "Mix Low Result B")), map_signif_level = TRUE, 
    test = "t.test", textsize = 6) + geom_signif(comparisons = list(c("Mix Low Result C", 
    "Mix Low Result D")), map_signif_level = TRUE, test = "t.test", textsize = 6) + 
    geom_signif(comparisons = list(c("Mix Low Result B", "Mix Low Result D")), map_signif_level = TRUE,
                test = "t.test", textsize = 6, y_position = 3) + 
geom_signif(comparisons = list(c("Mix High Result A", "Mix High Result B")), map_signif_level = TRUE, 
    test = "t.test", textsize = 6) + geom_signif(comparisons = list(c("Mix High Result C", 
    "Mix High Result D")), map_signif_level = TRUE, test = "t.test", textsize = 6) + 
    geom_signif(comparisons = list(c("Mix High Result B", "Mix High Result D")),
                map_signif_level = TRUE, test = "t.test", textsize = 6, y_position = 3) + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Plot m/z ranges integrated by XCMS
p.deltamz.box <- ggplot(data = summary, aes(x = experiment, y = deltamz)) + geom_boxplot(aes(group = experiment)) + 
    ggtitle("Mass ranges within feature groups in plasma samples") + xlab("Experiments") + 
    ylab("Log10 of mass ranges") + geom_signif(comparisons = list(c("Mix Low Result A", 
    "Mix Low Result B")), map_signif_level = TRUE, test = "t.test", textsize = 6) + 
    geom_signif(comparisons = list(c("Mix Low Result C", "Mix Low Result D")), map_signif_level = TRUE, 
        test = "t.test", textsize = 6) + geom_signif(comparisons = list(c("Mix Low Result B", 
    "Mix Low Result D")), map_signif_level = TRUE, test = "t.test", textsize = 6, 
    y_position = 0.5) + 
geom_signif(comparisons = list(c("Mix High Result A", "Mix High Result B")), map_signif_level = TRUE, 
    test = "t.test", textsize = 6) + geom_signif(comparisons = list(c("Mix High Result C", 
    "Mix High Result D")), map_signif_level = TRUE, test = "t.test", textsize = 6) + 
    geom_signif(comparisons = list(c("Mix High Result B", "Mix High Result D")), 
        map_signif_level = TRUE, test = "t.test", textsize = 6, y_position = 0.5) + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.deltamz.box <- ggplot_gtable(ggplot_build(p.deltamz.box))
p.deltamz.box$layout$clip[p.deltamz.box$layout$name == "panel"] <- "off"

## Save plots to harddrive
ggsave("cv.pdf", p.cv.box, width = 8, height = 6, device = "pdf")
ggsave("mz.pdf", p.deltamz.box, width = 8, height = 6, device = "pdf")
ggsave("rt.pdf", p.deltart.box, width = 8, height = 6, device = "pdf")
