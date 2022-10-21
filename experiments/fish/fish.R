####################################################################
# Your R session's cwd should be the root of the project
####################################################################
source("./src/plotutils.R")

####################################################################
# To create the dataset execute
####################################################################
# library(ForwardSearch)
# data('Fulton')
# save(Fulton, file = 'Fulton.Rda')
####################################################################
path <- './experiments/fish/'
load(paste0(path, 'Fulton.Rda'))
d <- Fulton
d$y <- d$q

mod1 <- d$Wed == 0
mod2 <- d$Stormy == 0
de_train_idx <- mod2 & mod1 # fair weather
do_train_idx <- !mod2 & mod1
d_test_idx <- !mod1

####################################################################
# visualize datasets
####################################################################

# visualize response
main_text_size <- 12
d_ <- d
d_$g <- "l"
d_[de_train_idx,]$g <- paste0(
  "In-Sample (Not Wednesday, Fair weather) (", sum(de_train_idx), " obs.)")
d_[do_train_idx,]$g <- paste0(
  "In-Sample (Not Wednesday, Stormy weather) (", sum(do_train_idx), " obs.)")
d_[d_test_idx,]$g <- paste0(
  "Out-of-Sample (Wednesday) (", sum(d_test_idx), " obs.)")
p <- ggplot(d_, aes(y, group = g, fill = g)) +
  geom_density(alpha = 0.5) +
  ylab('Approximate density') +
  xlab('log(quantity)') +
  scale_fill_discrete(name = "") +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(axis.text = element_text(size = main_text_size),
        legend.text = element_text(size = main_text_size))

# visualize covariate
p2 <- ggplot(d_, aes(p, group = g, fill = g)) +
  geom_density(alpha = 0.5) +
  ylab('Approximate density') +
  xlab('log(price)') +
  scale_fill_discrete(name = "") +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(axis.text = element_text(size = main_text_size),
        legend.text = element_text(size = main_text_size))

pp2 <- ggarrange(p2, p, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
save_plot(pp2, paste0(path, "fish_datasets.pdf"), dims = list(width = 15, height = 5))
####################################################################

source("./src/real.R")
source("./src/scr.R")

compute_estimator <- compute_scr
method <- (function(data) scr(data = data, length = 200))

set.seed(123)
result <- resample_test(d = d,
                        names = c("p"),
                        de_train_idx = de_train_idx,
                        do_train_idx = do_train_idx,
                        d_test_idx = d_test_idx,
                        method = method,
                        compute_estimator = compute_estimator,
                        resamples = 1000,
                        debug = TRUE,
                        folds = 3)

save(result, file = paste0(path, 'result.Rda'))

main_text_size <- 12
p1 <- ggplot(data.frame(x = result$cd - result$cr), aes(x = x)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
  geom_density(alpha = .2, fill = "#FF6666") +
  ylab('Approximate density') +
  xlab(expression("Risk"(hat(beta)[CD]) - "Risk"(hat(beta)[CR]))) +
  ggtitle('Risk difference in test set') +
  theme(axis.text = element_text(size = main_text_size), axis.title = element_text(size = main_text_size))
p2 <- result$cv$plot +
  ggtitle('Cross-validation in training set') +
  theme(legend.position = c(0.9, 0.9)) +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype = 'solid')) +
  theme(axis.text = element_text(size = main_text_size), axis.title = element_text(size = main_text_size))
q <- grid.arrange(p2, p1, ncol = 2)
save_plot(q, paste0(path, "fish_result.pdf"), dims = list(width = 10))







