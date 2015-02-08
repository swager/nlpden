
#' Density estimation by quadratic programming
#'
#' Estimates the density of a signal that has been
#' convolved with a known noise distribution.
#' For example, if X = mu + epsilon where epsilon is noise drawn from a known
#' noise_density_function, then nlpden will return an estimate of the
#' density f(X).
#'
#' @param noisy_signal the raw observations
#' @param noise_density_fn the density of the additive noise function
#' @param BIN_SZ The size of the histogram bins for the empirical distribution
#' @param MARGIN The amount of padding on each side of the data
#'
#' @export nlpden
#' @return density estimate of the form (x, f(x))
#'
#' @examples
#' mu = c(runif(N/2, min = -3, max = 3), rep(0, N/2))
#' X = mu + rnorm(N)
#' f.hat = nlpden(X)
#' plot(f.hat, type = "l")


nlpden = function (noisy_signal_raw,
    noise_density_fn = dnorm,
	BIN_SZ = 0.01,
    MARGIN = 1) {
  
  # Center the data first
  center = median(noisy_signal_raw);
  noisy_signal = noisy_signal_raw - center;  	
  
  # Preparation
  bins = seq(min(noisy_signal) - MARGIN, max(noisy_signal) + MARGIN, by = BIN_SZ);
  data.len = max(noisy_signal) - min(noisy_signal) + 2*MARGIN;
  plot.idx = bins[-1] + BIN_SZ/2;

  # Make histogram. The histogram is normalized to sum to 1.
  # In the last step, divide by the bin size to fix normalization
  data.hist = hist(noisy_signal, breaks=bins, plot=FALSE);
  data.raw = data.hist$counts/length(noisy_signal);
  bin.count = length(data.raw);
  idx = 0:(bin.count - 1);

  # Compute fft of noise function
  noise.density = sapply(plot.idx, noise_density_fn);
  noise.fft = fft(noise.density/sum(noise.density));

  # Do fft for data, and get l2-closest deconvolution
  data.fft = fft(data.raw);
  noise_limit = min(floor(bin.count/2), which(Mod(data.fft) > 100*Mod(noise.fft))) - 1;
  raw_partial_deconvolution.fft = data.fft[1:noise_limit]/noise.fft[1:noise_limit];
  noise_weight = Mod(noise.fft[1:noise_limit]);
  feasible_soln.fft = make_feasible(raw_partial_deconvolution.fft, noise_weight);
  
  # Pad central area with zeros, and get deconvolution  
  full_deconvolution.fft = rep(0, bin.count);
  full_deconvolution.fft[1:noise_limit] = feasible_soln.fft[1:noise_limit];
  full_deconvolution.fft[(bin.count - noise_limit + 2):bin.count] = feasible_soln.fft[(noise_limit + 1):(2*noise_limit - 1)];
  
  # Get density estimate
  density_estimate.fft = full_deconvolution.fft * noise.fft;
  density_estimate.unnormalized = Re(fft(density_estimate.fft, inverse = TRUE))/length(density_estimate.fft);
  density_estimate.normalized = density_estimate.unnormalized/BIN_SZ;
  
  # Tidy up
  ret = data.frame(
  	"x_axis" = plot.idx + center, # undo the original centering
  	"density_estimate" = density_estimate.normalized
  );
  
  return(ret);
}

# Helper function: Essentially just a wrapper around the solve.QP call.
#
# arr is an array of (complex, positive index) Fourier coefficients a_0, a_1, ..., a_n.
# weights indicates the penalty for changing any of them

make_feasible = function(arr, weight) {
	n = length(arr) - 1;
	a0 = 1;
	a = Re(arr[-1]);
	b = Im(arr[-1]); 
	x0_vec = c(a0, a, b);
	weight_vec = c(weight, weight[-1])^2;
	
	fourier_rep = matrix(0, 2*n + 1, 2*n + 1);
	fourier_rep[1, n+1] = 1;
	for(iter in 1:n) {
		fourier_rep[1 + iter, c(n + 1 - iter, n + 1 + iter)] = c(1, 1);
		fourier_rep[n + 1 + iter, c(n + 1 - iter, n + 1 + iter)] = complex(imaginary = 1) * c(-1, 1);
	}
	
	fft_rep = matrix(0, 2*n + 1, 2*n + 1);
	fft_rep[,1:(n + 1)] = fourier_rep[,(n + 1):(2*n + 1)];
	fft_rep[,(n + 2):(2*n + 1)] = fourier_rep[,1:n];
	
	basis_rep = matrix(0, 2*n + 1, 2*n + 1);
	for (iter in 1:(2*n + 1)) {
		basis_rep[iter,] = fft(fft_rep[iter,], inverse = TRUE)/(2*n + 1);
	}
	basis_rep = Re(basis_rep); # rest is just numeric error;
	
	qp.soln = quadprog::solve.QP(Dmat = diag(weight_vec[-1]), dvec = (x0_vec*weight_vec)[-1], Amat = basis_rep[-1,], bvec = -basis_rep[1,]);
	constrained_soln = c(1, qp.soln$solution);
	constrained_distn = t(basis_rep) %*% constrained_soln;
	new_arr = fft(constrained_distn); #this should be equivalent to constrained solution, if everything were perfect
	return(new_arr);
}
