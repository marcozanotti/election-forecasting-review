
# Installs & Load ---------------------------------------------------------

install.packages(lattice)

library(lattice)



# Simulations -------------------------------------------------------------

# define function
f <- deriv(~ x^2 + y^2 + 1.5 * x * y, c("x", "y"), function.arg = TRUE)

# plot function
f <- deriv(~ x^2 + y^2 + 1.5 * x * y, c("x", "y"), function.arg = TRUE)
xpts <- seq(-3, 3, len = 40)
ypts <- seq(-3, 3, len = 40)
feval <- outer(xpts, ypts, f)
wireframe(feval, drape = TRUE, col.regions = rainbow(100), zlab = "feval")

# SD simulation
x0 <- -2.5                           
y0 <- 1.2
sd_res <- steepest_descent(x0, y0, f, tol = 0.01, max_iter = 100)
sd_res
plot_search(x0, y0, f, sd_res)

# CG simulation
cg_res <- conjugate_gradient(x0, y0, f, tol = 0.01, max_iter = 100)
cg_res
plot_search(x0, y0, f, cg_res)

# combine results
plot_search(x0, y0, f, list(cg_res, sd_res))



# Linear Regression -------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(gapminder)
library(bench)
source("R/utils.R")

gap_clean <- gapminder::gapminder %>% 
	tidyr::pivot_wider(names_from = continent, values_from = continent) %>% 
	dplyr::mutate(dplyr::across(Asia:Oceania, ~ ifelse(is.na(.), 0, 1))) %>% 
	dplyr::mutate(dplyr::across(pop:gdpPercap, ~ as.numeric(scale(.)))) %>% 
	dplyr::mutate(Intercept = 1, .before = pop) %>% 
	dplyr::select(-country, -year, -Africa)

gap_clean
y_vec <- as.matrix(gap_clean$lifeExp)
x_mat <- as.matrix(gap_clean[, -1])

# compute betas coefficients
mod_cg <- conjugate_gradient_lm(y_vec, x_mat)
mod_lm <- lm(lifeExp ~ . - 1, data = gap_clean)
results <- data.frame(
	"beta_CG" = mod_cg %>% round(5),
	"beta_lm" = coefficients(mod_lm) %>% round(5)
)
results
all.equal(results$beta_CG, results$beta_lm)

# computational results
comp_res <- bench::mark(
	lm(lifeExp ~ . - 1, data = gap_clean),
	conjugate_gradient_lm(y_vec, x_mat),
	iterations = 1000, 
	check = FALSE
) %>% 
	dplyr::mutate(expression = c("LM", "GC")) %>% 
	dplyr::select(expression, median, mem_alloc, n_itr) %>% 
	dplyr::rename(method = expression)
comp_res




