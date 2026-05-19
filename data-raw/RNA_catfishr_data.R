# Creating RNA Space Test Data

RNA_Ref <- cbind(
  PC1 = rnorm(200, mean = 0.7, sd = 1),
  PC2 = rnorm(200, mean = 0.7, sd = 1),
  PC3 = rnorm(200, mean = 0.7, sd = 1)
)

RNA_Query_A <- cbind(
  PC1 = rnorm(600, mean = 3.4, sd = 2),
  PC2 = rnorm(600, mean = 3.4, sd = 2),
  PC3 = rnorm(600, mean = 3.4, sd = 2)
)

RNA_Query_B <- cbind(
  PC1 = rnorm(600, mean = 5, sd = 1.2),
  PC2 = rnorm(600, mean = 5, sd = 1.2),
  PC3 = rnorm(600, mean = 5, sd = 1.2)
)

pca_matrix <- rbind(RNA_Ref, RNA_Query_A, RNA_Query_B)

cell_ids <- paste0("cell_", seq_len(nrow(pca_matrix)))
rownames(pca_matrix) <- cell_ids

clusters <- c(
  rep("RNA_Ref", nrow(RNA_Ref)),
  rep("RNA_Query_A", nrow(RNA_Query_A)),
  rep("RNA_Query_B", nrow(RNA_Query_B))
)

names(clusters) <- cell_ids

RNA_catfishr_data <- list(
  pca_matrix = pca_matrix,
  clusters = clusters,
  ref_clusters = c("RNA_Ref")
)

usethis::use_data(RNA_catfishr_data, overwrite = TRUE)
