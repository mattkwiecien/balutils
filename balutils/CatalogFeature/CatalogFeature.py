from balutils import Catalog

# Essentially an interface for the decorator components

class CatalogFeature():
    """
    Defines the operations that decorators can implement.
    """
    def applyTo(self, catalog: Catalog) -> None:
        pass

