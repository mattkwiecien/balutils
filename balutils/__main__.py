from balutils.CatalogFactory import CatalogFactory
from balutils import Catalog


if __name__ == 'main':

    # pass in file name, list of properties that apply to this catalog.

    # fits catalog
    catalog = CatalogFactory.Create('myFileName', ['fits'])
    # fits, gold default, detection
    catalog = CatalogFactory.Create('myFileName', ['fits', 'gold_default', 'detection'])
    # etc.