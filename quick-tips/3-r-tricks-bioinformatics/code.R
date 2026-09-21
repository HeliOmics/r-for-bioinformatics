---
  title: "3R tricks in bioinformatics"
author: "Sanzida Akhter Anee"
date: "`r Sys.Date()`"
output: html_document
---
  
  
  # Step 1: Install & load packages
  

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma"))

library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)




# Step 2: Download the GEO dataset


gse <- getGEO("GSE19804", GSEMatrix = TRUE)
eset <- gse[[1]]   # ExpressionSet



# Step 3: Extract expression and metadata


expr <- exprs(eset)        # expression matrix (genes x samples)
meta <- pData(eset)        # sample metadata

dim(expr)
head(meta[, 1:5])


# Step 4: Explore date set


names(meta)
table(meta$group)
table(meta$`tissue:ch1`)





# Step 5: Define groups (condition labels)


# Example: using a tissue.ch column
# (Check names(meta) to confirm correct column)

# Check structure
str(meta)

# Create group properly
meta$group <- ifelse(
  meta$`tissue:ch1` == "lung cancer",
  "Tumor",
  "Normal"
)

meta$group <- factor(meta$group)

# Verify
table(meta$group)


# Step 6: Build design matrix and  git remote add origin git@github.com:sanzidaanee/youtube-shorts.git
git branch -M main
git push -u origin main
it linear model (limma)


library(limma)

design <- model.matrix(~ group, data = meta)

fit <- lmFit(expr, design)
fit <- eBayes(fit)

topTable(fit, coef = "groupTumor")




# Step 7:Extract differential expression results


results <- topTable(fit, coef = "groupTumor", number = Inf)

head(results)



#Step 8: Clean results table


genes <- results %>%
  mutate(gene = rownames(results)) %>%
  select(gene, logFC, P.Value, adj.P.Val)

head(genes)



#Step 9: Filter significant genes


sig_genes <- genes %>%
  filter(adj.P.Val < 0.05)

head(sig_genes)


# Step 10: Volcano plot (standard bioinformatics visualization)


ggplot(genes, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-log10(P-value)"
  )




# Step 11: Bar plot of top genes


ggplot(top_genes, aes(x = reorder(gene, logFC), y = logFC, fill = logFC)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  labs(
    title = "Top Differentially Expressed Genes",
    x = "Gene",
    y = "Log2 Fold Change"
  )






