# SD step n
steepest_descent_n <- function(x, y, f) {
	
	fval <- f(x, y) # evaluate f at x,y
	fg <- drop(attr(fval, "gradient")) # get the gradient
	p <- -fg # steepest descent direction
	
	# alpha function
	f_sub <- function(alpha) {
		ff <- f(x + alpha * p[1], y + alpha * p[2])
		return(as.numeric(ff))
	}
	alpha <- optimize(f_sub, c(0, 4))$minimum # minimize wrt alpha
	x1 <- x + alpha * p[1] # get x new
	y1 <- y + alpha * p[2] # get y new
	
	res <- list("x" = c(x1, y1),	"p" = p, "alpha" = alpha)
	return(res)
	
}

# SD full
steepest_descent <- function(x0, y0, f, tol = 0.01, max_iter = 100) {
	
	# initial step
	norm_n <- 1
	res_mat <- matrix(NA, max_iter, 5)
	n_iter <- 1
	res_n <- steepest_descent_n(x0, y0, f)
	res_mat[n_iter, ] <- unname(c(res_n$x, res_n$p, res_n$alpha))
	
	while (norm_n > tol) {
		res_n <- steepest_descent_n(res_n$x[1], res_n$x[2], f)
		n_iter <- n_iter + 1
		res_mat[n_iter, ] <- unname(c(res_n$x, res_n$p, res_n$alpha))
		norm_n <- norm(t(res_mat[(n_iter - 1):n_iter, 1:2]))
	}
	
	res_mat <- res_mat[!is.na(res_mat[,1]),]
	colnames(res_mat) <- c("x", "y", "px", "py", "alpha")
	return(res_mat)
	
}

# CG step n
conjugate_gradient_n <- function(x, y, p0, f) {
	
	fval <- f(x, y)
	fg <- drop(attr(fval, "gradient"))
	beta <- drop(crossprod(fg) / crossprod(p0)) # Fletcher-Reeves
	p <- -fg + beta * p0 # conjugate gradient direction
	
	# alpha function
	f_sub <- function(alpha) {
		ff <- f(x + alpha * p[1], y + alpha * p[2])
		return(as.numeric(ff))
	}
	alpha <- optimize(f_sub, c(0, 4))$minimum # minimize wrt alpha
	x1 <- x + alpha * p[1] # get x new
	y1 <- y + alpha * p[2] # get y new
	
	res <- list("x" = c(x1, y1),	"p" = p, "alpha" = alpha)
	return(res)
	
}

# CG full
conjugate_gradient <- function(x0, y0, f, tol = 0.01, max_iter = 100) {
	
	# initial step
	norm_n <- 1
	res_mat <- matrix(NA, max_iter, 5)
	n_iter <- 1
	res_n <- steepest_descent_n(x0, y0, f)
	res_mat[n_iter, ] <- unname(c(res_n$x, res_n$p, res_n$alpha))
	
	while (norm_n > tol) {
		res_n <- conjugate_gradient_n(res_n$x[1], res_n$x[2], res_n$p, f)
		n_iter <- n_iter + 1
		res_mat[n_iter, ] <- unname(c(res_n$x, res_n$p, res_n$alpha))
		norm_n <- norm(t(res_mat[(n_iter - 1):n_iter, 1:2]))
	}
	
	res_mat <- res_mat[!is.na(res_mat[,1]),]
	colnames(res_mat) <- c("x", "y", "px", "py", "alpha")
	return(res_mat)
	
}

# function to plot the line search
plot_search <- function(x0, y0, f, results) {
	
	x_start <- x0
	y_start <- y0
	n <- 40
	xpts <- seq(-3, 1, len = n)
	ypts <- seq(-1, 3, len = n)
	gr <- expand.grid(x = xpts, y = ypts)
	feval <- with(gr, f(x, y))
	z <- matrix(feval, nrow = n, ncol = n)
	
	par(mar = c(5, 4, 1, 1))
	contour(xpts, ypts, z, nlevels = 20)
	points(x0, y0, pch = 19, cex = 2)
	
	if (!is.list(results)) {
		xopt <- results[, 1]
		yopt <- results[, 2]
		for (i in 1:length(xopt)) {
			x1 <- xopt[i]
			y1 <- yopt[i]
			segments(x0, y0, x1, y1, lwd = 1)
			x0 <- x1
			y0 <- y1
		}
	} else {
		n_res <- length(results)
		for (nn in 1:n_res) {
			xopt <- results[[nn]][, 1]
			yopt <- results[[nn]][, 2]
			for (i in 1:length(xopt)) {
				x1 <- xopt[i]
				y1 <- yopt[i]
				segments(x0, y0, x1, y1, lwd = 2, lty = nn, col = nn)
				x0 <- x1
				y0 <- y1
			}	
			x0 <- x_start
			y0 <- y_start
		}
	}
	
	return(invisible(NULL))
	
}

# function to compute regression betas via CG optimization
conjugate_gradient_lm <- function(y_vec, x_mat) {
	
	beta_vec <- matrix(0, nrow = ncol(x_mat), ncol = 1) # initial guess for betas
	A <- t(x_mat) %*% x_mat # derive matrix A
	b <- t(x_mat) %*% y_vec # derive vector b
	fg <- A %*% beta_vec - b # residuals = gradient of f
	p <- -fg # p search direction
	k <- 0 # number of iterations
	
	while (norm(fg, "2") > 0.01) {
		
		alpha <- as.numeric((t(fg) %*% fg) / (t(p) %*% A %*% p)) # step length alpha
		beta_vec <- beta_vec + alpha * p # update beta coefficients
		fg1 <- fg + alpha * A %*% p # update gradient of f
		beta1 <- as.numeric((t(fg1) %*% fg1) / (t(fg) %*% fg)) # calculate beta for search direction
		p1 <- -fg1 + beta1 * p # update search direction (p)
		fg <- fg1
		p <- p1
		k <- k + 1
	}
	
	return(beta_vec)
	
}
