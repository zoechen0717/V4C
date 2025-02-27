# write your silhouette score unit tests here
import numpy as np
import pytest
from sklearn.cluster import KMeans as SklearnKMeans
from sklearn.metrics import silhouette_samples
from cluster import Silhouette

def test_silhouette():
    # Generate test data
    np.random.seed(42)
    X = np.random.rand(100, 2) * 10  # 100 points in 2D space
    k = 3

    # Fit sklearn KMeans for labels
    sklearn_kmeans = SklearnKMeans(n_clusters=k, random_state=42, n_init=10)
    sklearn_kmeans.fit(X)
    labels = sklearn_kmeans.labels_

    # Compute silhouette scores using sklearn
    sklearn_silhouette = silhouette_samples(X, labels)

    # Compute silhouette scores using custom implementation
    custom_silhouette = Silhouette().score(X, labels)

    # Check if silhouette scores are similar
    assert np.allclose(custom_silhouette, sklearn_silhouette, atol=1e-2)

    # Edge Cases
    with pytest.raises(ValueError):
        Silhouette().score(np.array([1, 2, 3]), labels)  # Input must be 2D
    with pytest.raises(ValueError):
        Silhouette().score(X, np.array([[0, 1], [1, 2]]))  # Labels must be 1D
    with pytest.raises(ValueError):
        Silhouette().score(X, np.array([0, 1, 2]))  # Labels length must match X.shape[0]
