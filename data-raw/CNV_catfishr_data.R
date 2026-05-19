# Creating CNV Space Test Data

CNV_Ref <- cbind(
  PC1 = rnorm(300, mean = 1, sd = 1),
  PC2 = rnorm(300, mean = 1, sd = 1),
  PC3 = rnorm(300, mean = 1, sd = 1)
)

CNV_Query_X <- cbind(
  PC1 = rnorm(500, mean = 1.5, sd = 1.2),
  PC2 = rnorm(500, mean = 1.5, sd = 1.2),
  PC3 = rnorm(500, mean = 1.5, sd = 1.2)
)

CNV_Query_Y <- cbind(
  PC1 = rnorm(400, mean = 10, sd = 1),
  PC2 = rnorm(400, mean = 10, sd = 1),
  PC3 = rnorm(400, mean = 10, sd = 1)
)

CNV_Query_Z <- cbind(
  PC1 = rnorm(200, mean = 7, sd = 1),
  PC2 = rnorm(200, mean = 7, sd = 1),
  PC3 = rnorm(200, mean = 7, sd = 1)
)

pca_matrix <- rbind(CNV_Ref, CNV_Query_X, CNV_Query_Y, CNV_Query_Z)

cell_ids <- paste0("cell_", seq_len(nrow(pca_matrix)))
rownames(pca_matrix) <- cell_ids

clusters <- c(
  rep("CNV_Ref", nrow(CNV_Ref)),
  rep("CNV_Query_X", nrow(CNV_Query_X)),
  rep("CNV_Query_Y", nrow(CNV_Query_Y)),
  rep("CNV_Query_Z", nrow(CNV_Query_Z))
)

names(clusters) <- cell_ids

CNV_catfishr_data <- list(
  pca_matrix = pca_matrix,
  clusters = clusters,
  ref_clusters = c("CNV_Ref")
)

usethis::use_data(CNV_catfishr_data, overwrite = TRUE)
