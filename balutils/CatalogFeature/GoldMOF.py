from balutils.CatalogFeature import SimpleCatalog
from balutils import Catalog


# could definitely abstract one "gold" catalog that defers to subclass for the column name
# specific to SOF... not sure if its worth the added complexity.

class Gold(SimpleCatalog):
    """
    Adds Gold catalog functionality to the catalog.
    """

    _gold_cut_cols = [
        'flags_foreground', 'flags_badregions', 'flags_footprint', 'meas_FLAGS_GOLD_MOF_ONLY'
    ]

    def applyTo(self, catalog: Catalog) -> None:
        catalog._check_for_cols(_gold_cut_cols)

        gold_cuts = np.where(
            (catalog._cat['flags_foreground'] == 0) &
            (catalog._cat['flags_badregions'] < 2) &
            (catalog._cat['flags_footprint'] == 1) &
            (catalog._cat['meas_FLAGS_GOLD_MOF_ONLY'] < 2)
        )

        catalog.apply_cut(gold_cuts)

        return
