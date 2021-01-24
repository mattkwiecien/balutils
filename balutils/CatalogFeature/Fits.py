from balutils.CatalogFeature import *

class Fits(DecoratorBase):
    """
    Adds fits catalog functionality to the catalog.
    """

    def operation(self) -> str:

        return self.component.operation()
