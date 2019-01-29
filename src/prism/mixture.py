import numpy as np

import cleanlog

from scipy.special import gammaln as logG
from scipy.special import digamma
from scipy.special import polygamma
from sklearn.mixture import GaussianMixture

logger = cleanlog.ColoredLogger(name='PRISM')

def trigamma(x):
    return polygamma(1, x)


class BetaBinomialMixture():
    """Beta-binomial mixture model for PRISM core algorithm.
    """
    def __init__(self, n_components=1, max_iter=10000, tol=1e-3, seed=None, verbose=False):
        self.n_components = n_components
        self.max_iter = max_iter
        self.tol = tol
        self.seed = seed
        self.verbose = verbose

    def _bebin_mle(self, ns, ks, ws, n_iter=4):
        """Given depths(ns), and counts(ks) find the most likely alpha and beta,
        which are the two parameters of beta-binomial distribution."""
        a, b = 1, 1
        N = ws.sum()
        p1_bar = np.exp(1 / N * np.sum(ws * (digamma(a + ks) - digamma(a + b + ns))))
        p2_bar = np.exp(1 / N * np.sum(ws * (digamma(b + (ns - ks)) - digamma(a + b + ns))))

        # Find good initial esimates.
        a = 1 / 2 * (1 - p2_bar) / (1 - p1_bar - p2_bar)
        b = 1 / 2 * (1 - p1_bar) / (1 - p1_bar - p2_bar)

        for _ in range(n_iter):
            p1_bar = np.exp(1 / N * np.sum(ws * (digamma(a + ks) - digamma(a + b + ns))))
            p2_bar = np.exp(1 / N * np.sum(ws * (digamma(b + (ns - ks)) - digamma(a + b + ns))))

            n_log_p1 = np.sum(ws * (digamma(a + ks) - digamma(a + b + ns)))
            n_log_p2 = np.sum(ws * (digamma(b + (ns - ks)) - digamma(a + b + ns)))

            q_1 = -N * trigamma(a)
            q_2 = -N * trigamma(b)
            g_1 = N * digamma(a + b) - N * digamma(a) + n_log_p1
            g_2 = N * digamma(a + b) - N * digamma(b) + n_log_p2
            z = N * trigamma(a + b)
            b_ = (g_1 / q_1 + g_2 / q_2) / (1 / z + 1 / q_1 + 1 / q_2)

            a = a - (g_1 - b_) / q_1
            b = b - (g_2 - b_) / q_2

        return a, b

    def _bebin_likelihood(self, n, k, a, b):
        return np.exp(self._bebin_loglikelihood(n, k, a, b))

    def _bebin_loglikelihood(self, n, k, a, b):
        tmp = logG(n + 1) + logG(k + a) + logG(n - k + b) + logG(a + b) - \
            (logG(k + 1) + logG(n - k + 1) + logG(a) + logG(b) + logG(n + a + b))
        return tmp.sum(axis=1)

    def _gmm_initialize(self, n, k):
        """Initialize alphas and betas by fitting gaussian mixture model roughly.
        Alphas and betas are correspondingly computed from means and variances of each component.
        """
        # Ratio of methylated reads.
        r = k / n
        r = r.reshape(self.n_data, self.n_dim)

        # Fit gaussian mixture model. It assumes that each dimension has its own variance.
        model = GaussianMixture(n_components=self.n_components, covariance_type='diag')
        model.fit(r)

        alphas, betas = [], []
        for i in range(self.n_components):
            # Compute alpha and beta from mean and variance of GMM fit.
            mu = model.means_[i]
            var = model.covariances_[i]

            alpha = ((1 - mu) / var - 1 / mu) * mu**2
            beta = alpha * (1 / mu - 1)

            alphas.append(alpha)
            betas.append(beta)
        return np.array(alphas), np.array(betas)

    def fit(self, n, k, headers):
        if self.seed:
            np.random.seed(self.seed)
        
        self.n_data, self.n_dim = len(n), 1 if n.ndim == 1 else len(n[0])

        # If 1d data are given, reshape them into 2d matrix.
        if self.n_dim == 1:
            n = n.reshape(self.n_data, self.n_dim)
            k = k.reshape(self.n_data, self.n_dim)

        self.depths, self.counts, self.headers = n, k, headers
        self.alphas_, self.betas_ = self._gmm_initialize(n, k)
        self.pi_ = np.ones(self.n_components) / self.n_components

        prev_weighted_loglikelihood, self.converged_ = np.inf, False

        for iteration in range(self.max_iter):
            # Compute likelihood of alpha and beta w.r.t. each data point.
            l = np.array([self._bebin_likelihood(n, k, a, b) for a, b in zip(self.alphas_, self.betas_)])

            # Compute posterior.
            weighted_likelihood = l * self.pi_.reshape([self.n_components, 1])
            curr_weighted_loglikelihood = np.log(weighted_likelihood.sum(axis=0)).sum()

            # Check for convergence.
            if np.abs(curr_weighted_loglikelihood - prev_weighted_loglikelihood) < self.tol:
                if self.verbose:
                    logger.debug('Met convergence criterion at iteration %d. Terminating.' % (iteration))
                self.converged_ = True
                break
            prev_weighted_loglikelihood = curr_weighted_loglikelihood

            w = weighted_likelihood / weighted_likelihood.sum(axis=0)

            # Compute relative size of clusters.
            self.pi_ = w.sum(axis=1) / self.n_data

            # Compute maximum likelihood estimates of alpha and beta.
            alphas, betas = [], []
            for i in range(self.n_components):
                alpha_list, beta_list = [], []
                for j in range(self.n_dim):
                    tmp_alpha, tmp_beta = self._bebin_mle(n[:, j], k[:, j], w[i])
                    alpha_list.append(tmp_alpha)
                    beta_list.append(tmp_beta)

                alphas.append(alpha_list)
                betas.append(beta_list)

            self.alphas_ = np.array(alphas)
            self.betas_ = np.array(betas)
        # END for iteration

        if not self.converged_:
            if self.verbose:
                logger.warning('The EM algorithm did not converge. Increase n_iter enough to ensure convergence.')

        self.log_likelihood_ = curr_weighted_loglikelihood

    def predict_proba(self, n, k):
        """Returns the posterior probabilities of each fingerprint epiloci for each cluster.

        :param list n: Depths of fingerprint epiloci.
        :param list k: Fingerprint pattern counts of fingerprint epiloci.

        :returns: Posterior probabilities of each fingerprint epiloci.
        """
        l = np.array([self._bebin_likelihood(n, k, a, b) for a, b in zip(self.alphas_, self.betas_)])
        weighted_likelihood = l * self.pi_.reshape([self.n_components, 1])
        return weighted_likelihood / weighted_likelihood.sum(axis=0)

    def _n_parameters(self):
        """Returns the number of parameters estimated while fitting the model.

        :returns: Number of parameters in the model.
        """
        return int(2 * self.n_dim * self.n_components + (self.n_components - 1))

    def get_weights(self):
        """Returns the list of cluster weights.
        Note that a cluster weight is computed as a sum of posterior probabilities that
        each of the data point will be assiged to that cluster.

        :returns: Cluster weights.
        """
        return self.pi_

    @property
    def means_(self):
        return np.array([a / (a + b) for a, b in zip(self.alphas_, self.betas_)])

    def get_means(self):
        """
        :returns: Cluster means.
        """
        return self.means_
    
    @property
    def dispersions_(self):
        return np.array([1 / (a + b + 1) for a, b in zip(self.alphas_, self.betas_)])

    def get_dispersions(self):
        """
        :returns: Cluster dispersions.
        """
        return self.dispersions_
    
    def bic(self):
        """
        :returns: Bayesian Information Criterion (BIC) value of the model.
        """
        return -2 * self.log_likelihood_ + np.log(self.n_data) * self._n_parameters()

    def get_n_dimensions(self):
        """
        :returns: Number of dimensions.
        """
        return self.n_dim

    def get_n_components(self):
        """
        :returns: Number of clusters.
        """
        return self.n_components

    def get_depths(self):
        """
        :returns: Depths of fingerprint epiloci used for fitting the model.
        """
        return self.depths
    
    def get_counts(self):
        """
        :returns: Fingerprint pattern counts of fingerprint epiloci used for fitting the model.
        """
        return self.counts
    
    def get_headers(self):
        """
        :returns: Headers of fingerprint epiloci used for fitting the model.
        """
        return self.headers
