set.seed(123)
library(data.table)

####################################################################
# Your R session's cwd should be the root of the project
####################################################################
source("./src/data.R") # load data generation methods
source("./src/datautil.R") # load data utilities
source("./src/measurements.R") # load metrics
source("./src/plotutils.R")
source("./src/selection.R") # load model selection methods
source("./src/scr.R")
path <- './experiments/cr_vs_cd/'

compare_shift <- function(target,
                          B,
                          shift,
                          n,
                          folds,
                          oshift = 100,
                          on = 10000) {

  # In-sample fit

  data <- data_gen(target = target, B = B, n = n, sd = shift)

  # Causal regularization

  estimands <- cross_validation(folds = folds, data = data, estimators = c(scr))

  l <- estimands[[1]]$opt
  m <- moments(data)

  cr <- compute_scr(m = m, l = l)
  cd <- compute_scr(m = m, l = 0)

  # Out-of-sample generation  
  odata <- gen(target = target, B = B, n = on, sd = oshift)

  return(list(

    # In-sample performance
    cr_in = in_sample_risk(data, beta = cr) / shift,
    cd_in = in_sample_risk(data, beta = cd) / shift,

    # Out-of-sample performance
    cr_out = risk(odata$X, odata$y, beta = cr) / oshift,
    cd_out = risk(odata$X, odata$y, beta = cd) / oshift

  ))

}

homogenize <- function(x, select, algo, shift) {
  s <- subset(x, select = select)
  s$algo <- algo
  s$shift <- shift
  names(s) <- c('insample', 'outsample', 'algo', 'shift')
  return(s)
}

collapse <- function(x, shift) {

  cr <- homogenize(x, c('cr_in', 'cr_out'), 'cr', shift)
  cd <- homogenize(x, c('cd_in', 'cd_out'), 'cd', shift)

  return(rbind(cr, cd))
}

compare_shifts <- function(target,
                           B,
                           n,
                           trials,
                           folds = 3,
                           shifts = c(0.5, 1, 1.5)) {


  trials_per_shift <- vector("list", length(shifts))

  for (k in 1:length(shifts)) {

    shift <- shifts[[k]]
    trials_ <- vector("list", trials)

    for (trial in 1:trials) {

      trials_[[trial]] <- compare_shift(target = target,
                                        B = B,
                                        shift = shift,
                                        folds = folds,
                                        n = n)

    }

    trials_per_shift[[k]] <- rbindlist(trials_, use.names = TRUE)

  }

  collapsed <- mapply(collapse, trials_per_shift, shifts, SIMPLIFY = FALSE)
  collapsed <- rbindlist(collapsed, use.names = TRUE)
  collapsed <- data.frame(collapsed)

  save(collapsed, file = paste(paste0(path, 'cr_vs_cd_'), n, '.Rda', sep = ''))

  return(collapsed)

}

target <- 4
B <- matrix(c(0, 0, 0, 0, 0, 0, 0,
              1, 0, 0, 0, 0, 0, 0,
              1, 1, 0, 0, 0, 0, 0,
              0, 1, 1, 0, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 0, 1, 0), nrow = 7, byrow = TRUE)

data_250 <- compare_shifts(target = target, B = B, n = 250, trials = 1000)

data_1000 <- compare_shifts(target, B, n = 1000, trials = 1000)

source("./src/plotutils.R")

data_250$n <- 250
data_1000$n <- 1000
df <- rbind(data_250, data_1000)
data_ <- df
data_$algo <- as.factor(data_$algo)
data_$shift <- as.factor(data_$shift)
data_$n <- as.factor(data_$n)
levels(data_$n) <- paste("Sample size =", levels(data_$n), sep = " ")
data_ <- data_ %>% dplyr::arrange(algo)
means <- aggregate(. ~ algo + shift + n, data_, median)
main_text_size <- 12
scientific10 <- function(x) parse(text = paste0("10^", as.character(log10(x))))
p <- ggplot() +
  geom_point(data_,
             mapping =
               aes(x = insample,
                   y = outsample,
                   shape = shift,
                   fill = algo,
                   color = algo),
             alpha = 0.15) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5)),
         shape = guide_legend(override.aes = list(alpha = 1, size = 5))) +
  geom_point(means,
             mapping =
               aes(x = insample,
                   y = outsample,
                   shape = shift,
                   fill = algo),
             alpha = 1,
             size = 5,
             show.legend = F) +
  facet_wrap(. ~ n, ncol = 2, scales = 'free') +
  scale_shape_manual(values = c(21, 24, 22)) +
  labs(color = "Algorithm", fill = "Algorithm", shape = "Initial Shift") +
  scale_fill_discrete(labels = c("causal Dantzig", "causal regularization")) +
  scale_color_discrete(labels = c("causal Dantzig", "causal regularization")) +
  xlab('Normalized in-sample risk (log scale)') +
  ylab('Normalized out-of-sample risk (log scale)') +
  scale_y_continuous(trans = 'log10', labels = scientific10, limits = c(NA, 1e2)) +
  scale_x_continuous(trans = 'log10', labels = scientific10, limits = c(NA, 1e3)) +
  theme(legend.position = "top") +
  theme(strip.text.x = element_text(size = main_text_size),
        legend.text = element_text(size = main_text_size),
        legend.title = element_text(size = main_text_size),
        axis.text = element_text(size = main_text_size),
        axis.title = element_text(size = main_text_size))
save_plot(p, paste0(path, "comparison2.pdf"), dims = list(width = 10, height = 5))



