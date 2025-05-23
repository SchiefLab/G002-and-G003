"""
This type stub file was generated by pyright.
"""

"""Statistical transformations for visualization.

This module is currently private, but is being written to eventually form part
of the public API.

The classes should behave roughly in the style of scikit-learn.

- All data-independent parameters should be passed to the class constructor.
- Each class should implement a default transformation that is exposed through
  __call__. These are currently written for vector arguments, but I think
  consuming a whole `plot_data` DataFrame and return it with transformed
  variables would make more sense.
- Some class have data-dependent preprocessing that should be cached and used
  multiple times (think defining histogram bins off all data and then counting
  observations within each bin multiple times per data subsets). These currently
  have unique names, but it would be good to have a common name. Not quite
  `fit`, but something similar.
- Alternatively, the transform interface could take some information about grouping
  variables and do a groupby internally.
- Some classes should define alternate transforms that might make the most sense
  with a different function. For example, KDE usually evaluates the distribution
  on a regular grid, but it would be useful for it to transform at the actual
  datapoints. Then again, this could be controlled by a parameter at  the time of
  class instantiation.

"""

class KDE:
    """Univariate and bivariate kernel density estimator."""

    def __init__(self, *, bw_method=..., bw_adjust=..., gridsize=..., cut=..., clip=..., cumulative=...) -> None:
        """Initialize the estimator with its parameters.

        Parameters
        ----------
        bw_method : string, scalar, or callable, optional
            Method for determining the smoothing bandwidth to use; passed to
            :class:`scipy.stats.gaussian_kde`.
        bw_adjust : number, optional
            Factor that multiplicatively scales the value chosen using
            ``bw_method``. Increasing will make the curve smoother. See Notes.
        gridsize : int, optional
            Number of points on each dimension of the evaluation grid.
        cut : number, optional
            Factor, multiplied by the smoothing bandwidth, that determines how
            far the evaluation grid extends past the extreme datapoints. When
            set to 0, truncate the curve at the data limits.
        clip : pair of numbers or None, or a pair of such pairs
            Do not evaluate the density outside of these limits.
        cumulative : bool, optional
            If True, estimate a cumulative distribution function. Requires scipy.

        """
        ...
    def define_support(
        self, x1, x2=..., weights=..., cache=...
    ):  # -> NDArray[floating[Any]] | tuple[NDArray[floating[Any]], NDArray[floating[Any]]]:
        """Create the evaluation grid for a given data set."""
        ...
    def __call__(
        self, x1, x2=..., weights=...
    ):  # -> tuple[NDArray[Unknown] | Unknown | Any, NDArray[floating[Any]] | tuple[NDArray[floating[Any]], NDArray[floating[Any]]]] | tuple[NDArray[float64] | Unknown | Any, NDArray[floating[Any]] | tuple[NDArray[floating[Any]], NDArray[floating[Any]]]]:
        """Fit and evaluate on univariate or bivariate data."""
        ...

class Histogram:
    """Univariate and bivariate histogram estimator."""

    def __init__(self, stat=..., bins=..., binwidth=..., binrange=..., discrete=..., cumulative=...) -> None:
        """Initialize the estimator with its parameters.

        Parameters
        ----------
        stat : str
            Aggregate statistic to compute in each bin.

            - `count`: show the number of observations in each bin
            - `frequency`: show the number of observations divided by the bin width
            - `probability` or `proportion`: normalize such that bar heights sum to 1
            - `percent`: normalize such that bar heights sum to 100
            - `density`: normalize such that the total area of the histogram equals 1

        bins : str, number, vector, or a pair of such values
            Generic bin parameter that can be the name of a reference rule,
            the number of bins, or the breaks of the bins.
            Passed to :func:`numpy.histogram_bin_edges`.
        binwidth : number or pair of numbers
            Width of each bin, overrides ``bins`` but can be used with
            ``binrange``.
        binrange : pair of numbers or a pair of pairs
            Lowest and highest value for bin edges; can be used either
            with ``bins`` or ``binwidth``. Defaults to data extremes.
        discrete : bool or pair of bools
            If True, set ``binwidth`` and ``binrange`` such that bin
            edges cover integer values in the dataset.
        cumulative : bool
            If True, return the cumulative statistic.

        """
        ...
    def define_bin_params(
        self, x1, x2=..., weights=..., cache=...
    ):  # -> dict[str, int | tuple[Any, Any]] | dict[str, NDArray[signedinteger[Any]] | NDArray[Any]] | dict[str, tuple[Unknown, ...]]:
        """Given data, return numpy.histogram parameters to define bins."""
        ...
    def __call__(
        self, x1, x2=..., weights=...
    ):  # -> tuple[Any | ndarray[Unknown, Unknown] | NDArray[float64] | NDArray[Any], NDArray[Any]] | tuple[Any | ndarray[Unknown, Unknown] | NDArray[float64], list[NDArray[floating[Any]]]]:
        """Count the occurrences in each bin, maybe normalize."""
        ...

class ECDF:
    """Univariate empirical cumulative distribution estimator."""

    def __init__(self, stat=..., complementary=...) -> None:
        """Initialize the class with its parameters

        Parameters
        ----------
        stat : {{"proportion", "count"}}
            Distribution statistic to compute.
        complementary : bool
            If True, use the complementary CDF (1 - CDF)

        """
        ...
    def __call__(self, x1, x2=..., weights=...):  # -> tuple[Any, Any]:
        """Return proportion or count of observations below each sorted datapoint."""
        ...

class EstimateAggregator:
    def __init__(self, estimator, errorbar=..., **boot_kws) -> None:
        """
        Data aggregator that produces an estimate and error bar interval.

        Parameters
        ----------
        estimator : callable or string
            Function (or method name) that maps a vector to a scalar.
        errorbar : string, (string, number) tuple, or callable
            Name of errorbar method (either "ci", "pi", "se", or "sd"), or a tuple
            with a method name and a level parameter, or a function that maps from a
            vector to a (min, max) interval.
        boot_kws
            Additional keywords are passed to bootstrap when error_method is "ci".

        """
        ...
    def __call__(self, data, var):  # -> Series[Unknown]:
        """Aggregate over `var` column of `data` with estimate and error interval."""
        ...
