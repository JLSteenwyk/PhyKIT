"""
Geological timescale data (ICS 2024).

Provides era, period, and epoch boundaries for plotting chronograms
with geological time bands.
"""

# Muted pastel color palette for geological epochs
# Inspired by ICS colors but softened for publication figures
EPOCH_COLORS = {
    "Holocene": "#FEF2CC",
    "Pleistocene": "#FFF2AE",
    "Pliocene": "#FFFF99",
    "Miocene": "#FFFF00",
    "Oligocene": "#FDC07A",
    "Eocene": "#FDB46C",
    "Paleocene": "#FDA75F",
    "Late Cretaceous": "#A6D84A",
    "Early Cretaceous": "#8CCD57",
    "Late Jurassic": "#B3E1EB",
    "Middle Jurassic": "#80CFD8",
    "Early Jurassic": "#42B2D0",
    "Late Triassic": "#BD8CC2",
    "Middle Triassic": "#B07FB7",
    "Early Triassic": "#A372AB",
}

PERIOD_COLORS = {
    "Quaternary": "#F9F97F",
    "Neogene": "#FFE619",
    "Paleogene": "#FD9A52",
    "Cretaceous": "#7FC64E",
    "Jurassic": "#34B2C9",
    "Triassic": "#812B92",
    "Permian": "#F04028",
    "Carboniferous": "#67A599",
    "Devonian": "#CB8C37",
    "Silurian": "#B3E1B6",
    "Ordovician": "#009270",
    "Cambrian": "#7FA056",
}

ERA_COLORS = {
    "Cenozoic": "#F2F91D",
    "Mesozoic": "#67C5CA",
    "Paleozoic": "#99C08D",
}

EPOCHS = [
    ("Holocene", 0.0117, 0),
    ("Pleistocene", 2.58, 0.0117),
    ("Pliocene", 5.333, 2.58),
    ("Miocene", 23.03, 5.333),
    ("Oligocene", 33.9, 23.03),
    ("Eocene", 56.0, 33.9),
    ("Paleocene", 66.0, 56.0),
    ("Late Cretaceous", 100.5, 66.0),
    ("Early Cretaceous", 145.0, 100.5),
    ("Late Jurassic", 163.5, 145.0),
    ("Middle Jurassic", 174.7, 163.5),
    ("Early Jurassic", 201.4, 174.7),
    ("Late Triassic", 237.0, 201.4),
    ("Middle Triassic", 247.2, 237.0),
    ("Early Triassic", 251.902, 247.2),
]

PERIODS = [
    ("Quaternary", 2.58, 0),
    ("Neogene", 23.03, 2.58),
    ("Paleogene", 66.0, 23.03),
    ("Cretaceous", 145.0, 66.0),
    ("Jurassic", 201.4, 145.0),
    ("Triassic", 251.902, 201.4),
    ("Permian", 298.9, 251.902),
    ("Carboniferous", 358.9, 298.9),
    ("Devonian", 419.2, 358.9),
    ("Silurian", 443.8, 419.2),
    ("Ordovician", 485.4, 443.8),
    ("Cambrian", 538.8, 485.4),
]

ERAS = [
    ("Cenozoic", 66.0, 0),
    ("Mesozoic", 251.902, 66.0),
    ("Paleozoic", 538.8, 251.902),
]


def get_timescale_for_range(root_age, level="auto"):
    """Return appropriate timescale intervals for a given root age.

    Parameters
    ----------
    root_age : float, Ma
    level : "epoch", "period", "era", or "auto"

    Returns
    -------
    list of (name, start_ma, end_ma), color_dict
    """
    if level == "auto":
        if root_age <= 66:
            level = "epoch"
        elif root_age <= 252:
            level = "period"
        else:
            level = "era"

    if level == "epoch":
        intervals = [(n, s, e) for n, s, e in EPOCHS if s > 0 and e < root_age * 1.1]
        # Include any epoch that overlaps the tree's range
        intervals = [(n, s, e) for n, s, e in EPOCHS if e < root_age * 1.05]
        return intervals, EPOCH_COLORS
    elif level == "period":
        intervals = [(n, s, e) for n, s, e in PERIODS if e < root_age * 1.05]
        return intervals, PERIOD_COLORS
    else:
        intervals = [(n, s, e) for n, s, e in ERAS if e < root_age * 1.05]
        return intervals, ERA_COLORS
