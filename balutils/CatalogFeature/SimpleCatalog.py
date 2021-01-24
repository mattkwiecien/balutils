from balutils.CatalogFeature import CatalogFeature
from balutils import Catalog


# I'm guessing that instead of "simple" the two base catalog types would be fits and H5 

class SimpleCatalog(CatalogFeature):
    """
    Base catalog features available.
    """
    def applyTo(self, catalog: Catalog) -> None:
        pass

