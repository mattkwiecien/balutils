from balutils.CatalogFeature import *
from balutils import Catalog
from astropy.table import Table
import fitsio

class Fits(DecoratorBase):
    """
    Adds fits catalog functionality to the catalog.
    """

    def applyTo(self, catalog: Catalog) -> None:
        
        # First, call the parent feature method
        self.parent.applyTo(catalog)

        # Then add my own logic
        fitsCatalog = Table(fitsio.read(Catalog.filename, columns=Catalog.cols))
        catalog._cat = fitsCatalog
        catalog.Nobjs = len(fitsCatalog)

        return
