from balutils.CatalogFeature import * 
from balutils import Catalog

class CatalogFactory():

    _featureLookup = {
        'fits': Fits,
        'gold_default': GoldDefault,
        'gold_SOF': GoldSOF,
        'gold_MOF': GoldMOF,
        'detection': Detection
    }

    def Create(self, filename, catalogTypes) -> Catalog:

        # this will probably end up being either an H5 or fits catalog, i think.
        catFeatures = SimpleCatalog()
        catalog = Catalog(filename)

        # loop through each type provided, and attach it's associated decorator
        for catalogType in catalogTypes:
            catFeatures = _featureLookup[catalogType](catFeatures)

        catFeatures.applyTo(catalog)

        return catalog
