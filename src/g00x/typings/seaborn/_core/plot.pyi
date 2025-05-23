"""
This type stub file was generated by pyright.
"""

import sys
from collections.abc import Callable, Generator
from contextlib import contextmanager
from typing import TYPE_CHECKING, Any, TypedDict

from matplotlib.axes import Axes
from matplotlib.figure import Figure, SubFigure
from seaborn._core.data import PlotData
from seaborn._core.moves import Move
from seaborn._core.scales import Scale
from seaborn._core.typing import (
    DataSource,
    Default,
    OrderSpec,
    VariableSpec,
    VariableSpecList,
)
from seaborn._marks.base import Mark
from seaborn._stats.base import Stat

"""The classes for specifying and compiling a declarative visualization."""
if TYPE_CHECKING: ...
if sys.version_info >= (3, 8): ...
else: ...
default = ...

class Layer(TypedDict, total=False):
    mark: Mark
    stat: Stat | None
    move: Move | list[Move] | None
    data: PlotData
    source: DataSource
    vars: dict[str, VariableSpec]
    orient: str
    legend: bool
    ...

class FacetSpec(TypedDict, total=False):
    variables: dict[str, VariableSpec]
    structure: dict[str, list[str]]
    wrap: int | None
    ...

class PairSpec(TypedDict, total=False):
    variables: dict[str, VariableSpec]
    structure: dict[str, list[str]]
    cross: bool
    wrap: int | None
    ...

@contextmanager
def theme_context(params: dict[str, Any]) -> Generator:
    """Temporarily modify specifc matplotlib rcParams."""
    ...

def build_plot_signature(cls):
    """
    Decorator function for giving Plot a useful signature.

    Currently this mostly saves us some duplicated typing, but we would
    like eventually to have a way of registering new semantic properties,
    at which point dynamic signature generation would become more important.

    """
    ...

@build_plot_signature
class Plot:
    """
    An interface for declaratively specifying statistical graphics.

    Plots are constructed by initializing this class and adding one or more
    layers, comprising a `Mark` and optional `Stat` or `Move`.  Additionally,
    faceting variables or variable pairings may be defined to divide the space
    into multiple subplots. The mappings from data values to visual properties
    can be parametrized using scales, although the plot will try to infer good
    defaults when scales are not explicitly defined.

    The constructor accepts a data source (a :class:`pandas.DataFrame` or
    dictionary with columnar values) and variable assignments. Variables can be
    passed as keys to the data source or directly as data vectors.  If multiple
    data-containing objects are provided, they will be index-aligned.

    The data source and variables defined in the constructor will be used for
    all layers in the plot, unless overridden or disabled when adding a layer.

    The following variables can be defined in the constructor:
        {known_properties}

    The `data`, `x`, and `y` variables can be passed as positional arguments or
    using keywords. Whether the first positional argument is interpreted as a
    data source or `x` variable depends on its type.

    The methods of this class return a copy of the instance; use chaining to
    build up a plot through multiple calls. Methods can be called in any order.

    Most methods only add information to the plot spec; no actual processing
    happens until the plot is shown or saved. It is also possible to compile
    the plot without rendering it to access the lower-level representation.

    """

    _data: PlotData
    _layers: list[Layer]
    _scales: dict[str, Scale]
    _shares: dict[str, bool | str]
    _limits: dict[str, tuple[Any, Any]]
    _labels: dict[str, str | Callable[[str], str]]
    _theme: dict[str, Any]
    _facet_spec: FacetSpec
    _pair_spec: PairSpec
    _figure_spec: dict[str, Any]
    _subplot_spec: dict[str, Any]
    _layout_spec: dict[str, Any]
    def __init__(self, *args: DataSource | VariableSpec, data: DataSource = ..., **variables: VariableSpec) -> None: ...
    def __add__(self, other): ...
    def on(self, target: Axes | SubFigure | Figure) -> Plot:
        """
        Provide existing Matplotlib figure or axes for drawing the plot.

        When using this method, you will also need to explicitly call a method that
        triggers compilation, such as :meth:`Plot.show` or :meth:`Plot.save`. If you
        want to postprocess using matplotlib, you'd need to call :meth:`Plot.plot`
        first to compile the plot without rendering it.

        Parameters
        ----------
        target : Axes, SubFigure, or Figure
            Matplotlib object to use. Passing :class:`matplotlib.axes.Axes` will add
            artists without otherwise modifying the figure. Otherwise, subplots will be
            created within the space of the given :class:`matplotlib.figure.Figure` or
            :class:`matplotlib.figure.SubFigure`.

        Examples
        --------
        .. include:: ../docstrings/objects.Plot.on.rst

        """
        ...
    def add(
        self,
        mark: Mark,
        *transforms: Stat | Mark,
        orient: str | None = ...,
        legend: bool = ...,
        data: DataSource = ...,
        **variables: VariableSpec
    ) -> Plot:
        """
        Specify a layer of the visualization in terms of mark and data transform(s).

        This is the main method for specifying how the data should be visualized.
        It can be called multiple times with different arguments to define
        a plot with multiple layers.

        Parameters
        ----------
        mark : :class:`Mark`
            The visual representation of the data to use in this layer.
        transforms : :class:`Stat` or :class:`Move`
            Objects representing transforms to be applied before plotting the data.
            Currently, at most one :class:`Stat` can be used, and it
            must be passed first. This constraint will be relaxed in the future.
        orient : "x", "y", "v", or "h"
            The orientation of the mark, which also affects how transforms are computed.
            Typically corresponds to the axis that defines groups for aggregation.
            The "v" (vertical) and "h" (horizontal) options are synonyms for "x" / "y",
            but may be more intuitive with some marks. When not provided, an
            orientation will be inferred from characteristics of the data and scales.
        legend : bool
            Option to suppress the mark/mappings for this layer from the legend.
        data : DataFrame or dict
            Data source to override the global source provided in the constructor.
        variables : data vectors or identifiers
            Additional layer-specific variables, including variables that will be
            passed directly to the transforms without scaling.

        Examples
        --------
        .. include:: ../docstrings/objects.Plot.add.rst

        """
        ...
    def pair(
        self, x: VariableSpecList = ..., y: VariableSpecList = ..., wrap: int | None = ..., cross: bool = ...
    ) -> Plot:
        """
        Produce subplots by pairing multiple `x` and/or `y` variables.

        Parameters
        ----------
        x, y : sequence(s) of data vectors or identifiers
            Variables that will define the grid of subplots.
        wrap : int
            When using only `x` or `y`, "wrap" subplots across a two-dimensional grid
            with this many columns (when using `x`) or rows (when using `y`).
        cross : bool
            When False, zip the `x` and `y` lists such that the first subplot gets the
            first pair, the second gets the second pair, etc. Otherwise, create a
            two-dimensional grid from the cartesian product of the lists.

        Examples
        --------
        .. include:: ../docstrings/objects.Plot.pair.rst

        """
        ...
    def facet(
        self,
        col: VariableSpec = ...,
        row: VariableSpec = ...,
        order: OrderSpec | dict[str, OrderSpec] = ...,
        wrap: int | None = ...,
    ) -> Plot:
        """
        Produce subplots with conditional subsets of the data.

        Parameters
        ----------
        col, row : data vectors or identifiers
            Variables used to define subsets along the columns and/or rows of the grid.
            Can be references to the global data source passed in the constructor.
        order : list of strings, or dict with dimensional keys
            Define the order of the faceting variables.
        wrap : int
            When using only `col` or `row`, wrap subplots across a two-dimensional
            grid with this many subplots on the faceting dimension.

        Examples
        --------
        .. include:: ../docstrings/objects.Plot.facet.rst

        """
        ...
    def scale(self, **scales: Scale) -> Plot:
        """
        Specify mappings from data units to visual properties.

        Keywords correspond to variables defined in the plot, including coordinate
        variables (`x`, `y`) and semantic variables (`color`, `pointsize`, etc.).

        A number of "magic" arguments are accepted, including:
            - The name of a transform (e.g., `"log"`, `"sqrt"`)
            - The name of a palette (e.g., `"viridis"`, `"muted"`)
            - A tuple of values, defining the output range (e.g. `(1, 5)`)
            - A dict, implying a :class:`Nominal` scale (e.g. `{"a": .2, "b": .5}`)
            - A list of values, implying a :class:`Nominal` scale (e.g. `["b", "r"]`)

        For more explicit control, pass a scale spec object such as :class:`Continuous`
        or :class:`Nominal`. Or use `None` to use an "identity" scale, which treats data
        values as literally encoding visual properties.

        Examples
        --------
        .. include:: ../docstrings/objects.Plot.scale.rst

        """
        ...
    def share(self, **shares: bool | str) -> Plot:
        """
        Control sharing of axis limits and ticks across subplots.

        Keywords correspond to variables defined in the plot, and values can be
        boolean (to share across all subplots), or one of "row" or "col" (to share
        more selectively across one dimension of a grid).

        Behavior for non-coordinate variables is currently undefined.

        Examples
        --------
        .. include:: ../docstrings/objects.Plot.share.rst

        """
        ...
    def limit(self, **limits: tuple[Any, Any]) -> Plot:
        """
        Control the range of visible data.

        Keywords correspond to variables defined in the plot, and values are a
        `(min, max)` tuple (where either can be `None` to leave unset).

        Limits apply only to the axis; data outside the visible range are
        still used for any stat transforms and added to the plot.

        Behavior for non-coordinate variables is currently undefined.

        Examples
        --------
        .. include:: ../docstrings/objects.Plot.limit.rst

        """
        ...
    def label(self, *, title=..., **variables: str | Callable[[str], str]) -> Plot:
        """
        Control the labels and titles for axes, legends, and subplots.

        Additional keywords correspond to variables defined in the plot.
        Values can be one of the following types:

        - string (used literally; pass "" to clear the default label)
        - function (called on the default label)

        For coordinate variables, the value sets the axis label.
        For semantic variables, the value sets the legend title.
        For faceting variables, `title=` modifies the subplot-specific label,
        while `col=` and/or `row=` add a label for the faceting variable.
        When using a single subplot, `title=` sets its title.

        Examples
        --------
        .. include:: ../docstrings/objects.Plot.label.rst


        """
        ...
    def layout(self, *, size: tuple[float, float] | Default = ..., engine: str | None | Default = ...) -> Plot:
        """
        Control the figure size and layout.

        .. note::

            Default figure sizes and the API for specifying the figure size are subject
            to change in future "experimental" releases of the objects API. The default
            layout engine may also change.

        Parameters
        ----------
        size : (width, height)
            Size of the resulting figure, in inches. Size is inclusive of legend when
            using pyplot, but not otherwise.
        engine : {{"tight", "constrained", None}}
            Name of method for automatically adjusting the layout to remove overlap.
            The default depends on whether :meth:`Plot.on` is used.

        Examples
        --------
        .. include:: ../docstrings/objects.Plot.layout.rst

        """
        ...
    def theme(self, *args: dict[str, Any]) -> Plot:
        """
        Control the default appearance of elements in the plot.

        .. note::

            The API for customizing plot appearance is not yet finalized.
            Currently, the only valid argument is a dict of matplotlib rc parameters.
            (This dict must be passed as a positional argument.)

            It is likely that this method will be enhanced in future releases.

        Matplotlib rc parameters are documented on the following page:
        https://matplotlib.org/stable/tutorials/introductory/customizing.html

        Examples
        --------
        .. include:: ../docstrings/objects.Plot.theme.rst

        """
        ...
    def save(self, loc, **kwargs) -> Plot:
        """
        Compile the plot and write it to a buffer or file on disk.

        Parameters
        ----------
        loc : str, path, or buffer
            Location on disk to save the figure, or a buffer to write into.
        kwargs
            Other keyword arguments are passed through to
            :meth:`matplotlib.figure.Figure.savefig`.

        """
        ...
    def show(self, **kwargs) -> None:
        """
        Compile the plot and display it by hooking into pyplot.

        Calling this method is not necessary to render a plot in notebook context,
        but it may be in other environments (e.g., in a terminal). After compiling the
        plot, it calls :func:`matplotlib.pyplot.show` (passing any keyword parameters).

        Unlike other :class:`Plot` methods, there is no return value. This should be
        the last method you call when specifying a plot.

        """
        ...
    def plot(self, pyplot: bool = ...) -> Plotter:
        """
        Compile the plot spec and return the Plotter object.
        """
        ...

class Plotter:
    """
    Engine for compiling a :class:`Plot` spec into a Matplotlib figure.

    This class is not intended to be instantiated directly by users.

    """

    _data: PlotData
    _layers: list[Layer]
    _figure: Figure
    def __init__(self, pyplot: bool, theme: dict[str, Any]) -> None: ...
    def save(self, loc, **kwargs) -> Plotter: ...
    def show(self, **kwargs) -> None:
        """
        Display the plot by hooking into pyplot.

        This method calls :func:`matplotlib.pyplot.show` with any keyword parameters.

        """
        ...
