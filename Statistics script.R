#PGx analysis AGORA-TM

#Hardy-Weinberg
install.packages("genetics")
library("genetics")


data <- c(rep("A/A", 9), rep("A/G", 12), rep("G/G", 21))
g <- genotype(data)
HWE.test(g)

#Fisher's exact test

star0_count_pop1 <- 200
total_count1 <- 500
star0_count_pop2 <- 450
total_count2 <- 500
star0_count_pop3 <- 336
total_count3 <- 500
star0_count_pop4 <- 250
total_count4 <- 500
star0_count_pop5 <- 367
total_count5 <- 500
star0_count_pop6 <- 407
total_count6 <- 500

df0 <- data.frame(
  "pop1" = c(star0_count_pop1, total_count1 - star0_count_pop1),
  "pop2" = c(star0_count_pop2, total_count2 - star0_count_pop2),
  "pop3" = c(star0_count_pop3, total_count3 - star0_count_pop3),
  "pop4" = c(star0_count_pop4, total_count4 - star0_count_pop4),
  "pop5" = c(star0_count_pop5, total_count5 - star0_count_pop5),
  "pop6" = c(star0_count_pop6, total_count6 - star0_count_pop6),
  row.names = c("starallele", "other"),
  stringsAsFactors = FALSE
)

colnames(df0) <- c("pop1", "pop2", "pop3", "pop4", "pop5","pop6")
df0


test1 <- fisher.test(subset(df0, select = c("pop1", "pop2")))
test1$p.value
test2 <- fisher.test(subset(df0, select = c("pop1", "pop3")))
test2$p.value
test3 <- fisher.test(subset(df0, select = c("pop1", "pop4")))
test3$p.value
test4 <- fisher.test(subset(df0, select = c("pop1", "pop5")))
test4$p.value
test5 <- fisher.test(subset(df0, select = c("pop1", "pop6")))
test5$p.value

#Benjamini-Yekutieli p-value adjustment
gene_data <- read.csv("path/to/file")
adjusted_p<-p.adjust(gene_data$Pvalue, method = 'BY', 40)
adjusted_gene <- cbind(gene_data, adjusted_p)


