import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KernelDensity
from scipy.stats import norm


"""
 This code is used to generate synthetic data using Monte Carlo Tree Search (MCTS)
 and Data Augmentation Markov chain Monte Carlo (DA-MCMC).
 If the data contains a combination of discreet and continuous variables, its better
 to use ctgan model of the Synthetic Data Vault project. 
"""

class MCTS:
    def __init__(self, data):
        """Initialize with continuous numerical data"""
        self.scaler = StandardScaler()
        self.scaled_data = self.scaler.fit_transform(data)
        
    def generate(self,n_clusters=5, n_samples=1000):
        """MCTS-inspired generation using cluster-guided sampling"""
        # Step 1: Cluster data
        kmeans = KMeans(n_clusters)
        clusters = kmeans.fit_predict(self.scaled_data)
        
        synthetic_samples = []
        
        # Step 2: Sample from clusters
        for _ in range(n_samples):
            # Random cluster selection
            cluster_id = np.random.choice(np.unique(clusters))
            cluster_data = self.scaled_data[clusters == cluster_id]
            
            # Generate sample from cluster statistics
            if len(cluster_data) > 0:
                sample = np.mean(cluster_data, axis=0) + np.random.normal(
                    scale=0.1, 
                    size=cluster_data.shape[1]
                )
                synthetic_samples.append(sample)
        
        # Step 3: Rescale to original distribution
        return self.scaler.inverse_transform(np.array(synthetic_samples))



class DAMCMC:
    def __init__(self, data, n_chains=5, kernel='gaussian'):
        """
        Initialize with continuous numerical data
        Args:
            data: 2D numpy array (n_samples, n_features)
            n_chains: Number of parallel chains
            kernel: as described in the scikit-learn library. 
                Options: 'gaussian', 'tophat', 'epanechnikov', 'exponential', 'linear', 'cosine'.
                Default: 'gaussian'
        """
        self.data = data
        self.n_chains = n_chains
        self.kde = KernelDensity(kernel=kernel)
        self.kde.fit(data)  # Initial KDE model for quick evaluation
        
    def _fast_log_prob(self, x):
        """First-stage evaluation using KDE"""
        return self.kde.score_samples(x.reshape(1, -1))[0]
    
    def _exact_log_prob(self, x):
        """Second-stage evaluation (can be replaced with more complex models)"""
        # Simple Gaussian mixture for demonstration
        return np.log(0.5 * norm.pdf(x, 0, 1).mean() + 
                      0.5 * norm.pdf(x, 5, 1).mean())
    
    def generate(self, n_samples=1000, burn_in=200):
        """DA-MCMC generation with two-stage acceptance"""
        # Initialize chains
        chains = np.random.permutation(self.data)[:self.n_chains]
        samples = []
        
        for _ in range(burn_in + n_samples):
            for i in range(self.n_chains):
                # Proposal generation
                proposal = chains[i] + np.random.normal(0, 0.1, size=chains[i].shape)
                
                # First stage (KDE evaluation)
                log_alpha1 = self._fast_log_prob(proposal) - self._fast_log_prob(chains[i])
                if np.log(np.random.rand()) < log_alpha1:
                    
                    # Second stage (exact evaluation)
                    log_alpha2 = self._exact_log_prob(proposal) - self._exact_log_prob(chains[i])
                    if np.log(np.random.rand()) < log_alpha2:
                        chains[i] = proposal
            
            if _ >= burn_in:
                samples.append(chains[np.random.randint(self.n_chains)].copy())
                
        return np.array(samples)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ks_2samp, wasserstein_distance
from sklearn.metrics import mutual_info_score
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

class SyntheticDataEvaluator:
    def __init__(self, original_df: pd.DataFrame, synthetic_df: pd.DataFrame):
        """
        Initialize the evaluator with original and synthetic dataframes.
        """
        self.original_df = original_df
        self.synthetic_df = synthetic_df

    def basic_stat_comparison(self):
        """
        Basic Statistical Comparison
        Prints descriptive statistics side-by-side for original and synthetic data.
        """
        print("\n=== Basic Statistics ===")
        stats_comparison = pd.concat([
            self.original_df.describe().add_prefix('original_'),
            self.synthetic_df.describe().add_prefix('synthetic_')
        ], axis=1)
        print(stats_comparison)

    def distribution_similarity(self):
        """
        Distribution Comparison Metrics
        Computes Kolmogorov-Smirnov test and Wasserstein distance for each feature.
        """
        print("\n=== Distribution Similarity ===")
        distribution_metrics = []

        for col in self.original_df.columns:
            ks_stat, ks_p = ks_2samp(self.original_df[col], self.synthetic_df[col])
            
            # Wasserstein Distance (Earth Mover's Distance)
            wasserstein = wasserstein_distance(self.original_df[col], self.synthetic_df[col])

            distribution_metrics.append({
                'Feature': col,
                'KS Statistic': ks_stat,
                'KS p-value': ks_p,
                'Wasserstein Distance': wasserstein
            })

        dist_df = pd.DataFrame(distribution_metrics)
        print(dist_df)
        return dist_df

    def correlation_preservation(self):
        """
        Correlation Preservation
        Plots heatmaps of correlation matrices for original and synthetic data.
        """
        print("\n=== Correlation Preservation ===")
        fig, ax = plt.subplots(1, 2, figsize=(15, 6))
        sns.heatmap(self.original_df.corr(), annot=True, fmt=".2f", ax=ax[0], cmap='coolwarm')
        ax[0].set_title("Original Data Correlation")
        sns.heatmap(self.synthetic_df.corr(), annot=True, fmt=".2f", ax=ax[1], cmap='coolwarm')
        ax[1].set_title("Synthetic Data Correlation")
        plt.show()

    def mutual_info_comparison(self):
        """
        Mutual Information Comparison
        Calculates mutual information between feature pairs for original and synthetic data.
        """
        print("\n=== Mutual Information ===")
        mi_comparison = []

        cols = self.original_df.columns
        for i, col1 in enumerate(cols):
            for col2 in cols[i+1:]:  # avoid duplicate pairs and self-pairs
                mi_original = mutual_info_score(self.original_df[col1], self.original_df[col2])
                mi_synthetic = mutual_info_score(self.synthetic_df[col1], self.synthetic_df[col2])
                mi_comparison.append({
                    'Feature Pair': f"{col1}-{col2}",
                    'Original MI': mi_original,
                    'Synthetic MI': mi_synthetic,
                    'Absolute Difference': abs(mi_original - mi_synthetic)
                })

        mi_df = pd.DataFrame(mi_comparison).sort_values('Absolute Difference', ascending=False)
        print(mi_df.head(10))
        return mi_df

    def nearest_neighbor_analysis(self):
        """
        Nearest Neighbor Distance Analysis
        Compares nearest neighbor distances within original and between synthetic and original data.
        """
        print("\n=== Nearest Neighbor Analysis ===")
        scaler = StandardScaler()
        original_scaled = scaler.fit_transform(self.original_df)
        synthetic_scaled = scaler.transform(self.synthetic_df)

        # Real-to-Real distances
        nn_real = NearestNeighbors(n_neighbors=2).fit(original_scaled)
        distances_real, _ = nn_real.kneighbors(original_scaled)
        real_mean_dist = np.mean(distances_real[:, 1])  # exclude self-distance

        # Synthetic-to-Real distances
        nn_synth = NearestNeighbors(n_neighbors=1).fit(original_scaled)
        distances_synth, _ = nn_synth.kneighbors(synthetic_scaled)
        synth_mean_dist = np.mean(distances_synth)

        print(f"Mean nearest neighbor distance (Real-Real): {real_mean_dist:.4f}")
        print(f"Mean nearest neighbor distance (Synth-Real): {synth_mean_dist:.4f}")

        return real_mean_dist, synth_mean_dist

    def distribution_visualisation(self):
        """
        Visualisation
        Plots KDE distributions of original and synthetic data for each feature.
        """
        print("\n=== Distribution Visualization ===")
        for col in self.original_df.columns:
            plt.figure(figsize=(10, 4))
            sns.kdeplot(self.original_df[col], label='Original', linewidth=2)
            sns.kdeplot(self.synthetic_df[col], label='Synthetic', linestyle='--')
            plt.title(f"Distribution Comparison: {col}")
            plt.legend()
            plt.show()

    def evaluate_all(self):
        """
        Run all evaluation steps sequentially.
        """
        self.basic_stat_comparison()
        self.distribution_similarity()
        self.correlation_preservation()
        self.mutual_info_comparison()
        self.nearest_neighbor_analysis()
        self.distribution_visualisation()
