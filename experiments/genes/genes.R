####################################################################
# Your R session's cwd should be the root of the project
####################################################################
source("./src/scr.R")
source("./src/real.R")
source("./src/plotutils.R")
path <- './experiments/genes/'
load(paste0(path, 'Kemmeren.RData'))

# Description
# data$obs    data without interventions
# data$int    data with interventions
# data$intpos indices of the intervened genes

intervened_idx <- as.vector(data$intpos)
n_observed <- dim(data$obs)[1]
n_intervened <- dim(data$int)[1] # number of intervections

# target selection
not_intervened <- !(1:6170 %in% intervened_idx)
observed_mean <- colSums(data$obs[, not_intervened]) / n_observed
intervened_mean <- colSums(data$int[, not_intervened]) / n_intervened
diff <- abs(observed_mean - intervened_mean)
or <- order(diff, decreasing = TRUE)
target_idx <- or[1] # most shifted non-intervened gene

# covariate selection
observed_mean <- colSums(data$obs[, intervened_idx]) / n_observed
intervened_mean <- colSums(data$int[, intervened_idx]) / n_intervened
diff <- abs(observed_mean - intervened_mean)
or <- order(diff, decreasing = TRUE)

# most shifted intervened genes
sampled_idx <- or[c(2, 3, 4, 5, 6, 7, 9, 11, 12, 13, 14, 15)] # for scr
sampled_observations <- intervened_idx %in% sampled_idx

yo <- as.matrix(data$obs[, target_idx])
Xo <- as.matrix(data$obs[, sampled_idx])

Xe <- as.matrix(data$int[, sampled_idx])
ye <- as.matrix(data$int[, target_idx])

ggplot(data.frame(y = c(yo, ye),
                  g = c(rep('yo', length(yo)),
                        rep('ye', length(ye))))) +
  geom_density(aes(y, group = g, fill = g), alpha = 0.5)


d <- list(Xe = Xe, ye = ye, Xo = Xo, yo = yo)

# check that the coefficient path has not singularities
s <- scr(data = center_train_data(d)$data, gamma = seq(0, 1, length.out = 300))
plot_path(gamma = s$gamma, beta = s$beta, names = data$genenames[sampled_idx])

###########################################
# choose method
###########################################

out_split <- 25 # outside split
in_split <- 5 # inside split
method <- (function(data) scr(data = data, length = 300))
compute_estimator <- compute_scr

###########################################

rcd <- vector(mode = "numeric", length = out_split)
rcr <- vector(mode = "numeric", length = out_split)

set.seed(123)
out_ge <- split_idx(out_split, n = n_intervened)
writeLines(paste("Average test dataset size", sum(out_ge == 1)))
writeLines(paste("Average training dataset size", sum(out_ge != 1)))

data_train <- vector("list", in_split)
data_val <- vector("list", in_split)
plots <- vector("list", out_split)
for (j in 1:out_split) {

  writeLines(paste("out_split", j))

  dcv <- subset_data(d, sele = (out_ge != j))
  centering <- center_train_data(dcv)
  dcv <- centering$data

  data_test <- subset_data(d, sele = (out_ge == j))
  data_test$ye <- data_test$ye - centering$mean

  in_ge <- split_idx(in_split, n = length(dcv$ye))
  if (j == 1) {
    writeLines(paste("Average inner-test dataset size", sum(in_ge == 1)))
    writeLines(paste("Average inner-train dataset size", sum(in_ge != 1)))
  }

  for (i in 1:in_split) {

    writeLines(paste("in_split", i))

    data_train[[i]] <- subset_data(dcv, sele = (in_ge != i))
    data_val[[i]] <- subset_data(dcv, sele = (in_ge == i))

  }

  sol <- fit_cross_validation(
    folds = in_split,
    f = method,
    data_train = data_train,
    data_val = data_val)

  l <- sol$opt
  m <- moments(dcv)

  # plot the selection and cv of one of the folds
  if (j == 25) {
    # plot covariate path
    s <- scr(dcv, gamma = seq(0, 1, length.out = 200))
    p2 <- plot_path(gamma = s$gamma, beta = s$beta, names = data$genenames[sampled_idx])
    p2 <- p2 + geom_vline(xintercept = l, color = "red", linetype = "dashed")
    ticks <- c(0, 0.25, 0.5, 0.75, 1, l)
    label_ticks <- expression(
      hat(beta)[CD], '0.25', '0.5', '0.75', hat(beta)[OLS], hat(beta)[CR])
    idx <- order(ticks)
    p2 <- p2 + scale_x_continuous(breaks = ticks[idx], labels = label_ticks[idx])

    p1 <- sol$plot
  }

  # causal regularization
  rcr[[j]] <- risk(X = data_test$Xe, y = data_test$ye, beta = compute_estimator(m, l = l))

  # causal dantzig
  rcd[[j]] <- risk(X = data_test$Xe, y = data_test$ye, beta = compute_estimator(m, l = 0))

}

result <- data.frame(cr = rcr, cd = rcd)

main_text_size <- 12
p1 <- p1 +
  theme(legend.position = c(0.85, 0.88)) +
  theme(legend.background =
          element_rect(colour = 'black', fill = 'white', linetype = 'solid')) +
  ggtitle('Cross-validation in training set') +
  theme(axis.text = element_text(size = main_text_size),
        legend.text = element_text(size = main_text_size),
        axis.title = element_text(size = main_text_size))
p2 <- p2 +
  guides(color = guide_legend(ncol = 4)) +
  theme(legend.position = c(0.65, 0.22)) +
  theme(legend.background =
          element_rect(colour = 'black', fill = 'white', linetype = 'solid')) +
  ggtitle('Coefficient path in training set') +
  theme(legend.title = element_blank()) +
  theme(axis.text = element_text(size = main_text_size),
        legend.text = element_text(size = 6),
        legend.spacing.y = unit(0.05, 'inches'),
        legend.spacing.x = unit(0.05, 'inches'),
        legend.margin = margin(t = 0, r = 0.05, b = 0, l = 0.05, unit = "inches"),
        axis.title = element_text(size = main_text_size))
p3 <- ggplot(data.frame(x = result$cd - result$cr), aes(x = x)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", bins =) +
  geom_density(alpha = .2, fill = "#FF6666", kernel = "epanechnikov") +
  ylab('Approximate density') +
  xlab(expression("Risk"(hat(beta)[CD]) - "Risk"(hat(beta)[CR]))) +
  ggtitle('Cross-validated risk difference in test set') +
  theme(axis.text = element_text(size = main_text_size), axis.title = element_text(size = main_text_size))
q <- grid.arrange(p1, p2, p3, layout_matrix = rbind(c(1, 2), c(1, 3)))
save_plot(q, filename = paste0(path, "genes_fold.pdf"), dims = list(width = 10, height = 6))








