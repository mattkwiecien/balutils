from balutils.CatalogFeature import SimpleCatalog

# TODO

class Matched(SimpleCatalog):
    """
    Adds Matched catalog functionality to the catalog.
    """
    def applyTo(self, catalog: Catalog) -> None:
        
        # First, call the parent feature method
        self.parent.applyTo(catalog)
        # TODO implement this decorators logic
        return
