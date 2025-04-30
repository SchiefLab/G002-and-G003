import pandas as pd


def top_n_percent(df: pd.DataFrame, col: str, percent: float = 0.9) -> pd.Series:
    """Qunatile of a column.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe
    col : str
        Column name; should be a numeric column
    percent : float, optional
        Percentile, by default 0.9

    Returns
    -------
    pd.Series
        Series with the quantile value
    """
    return pd.Series({col: df[col].quantile(percent, interpolation="midpoint")})


def top_n_percent_group(df: pd.DataFrame, group: list[str] | str, col: str, percent: float = 0.9) -> pd.DataFrame:
    """Top N percent from a group.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe
    group : list[str] | str
        Group column(s)
    metric : str
        DF column name; should be a numeric column
    percent : float, optional
        Percentile, by default 0.9

    Returns
    -------
    pd.DataFrame
        Dataframe with the quantile value
    """
    return df.groupby(group).apply(lambda x: top_n_percent(x, col=col, percent=percent)).reset_index()
