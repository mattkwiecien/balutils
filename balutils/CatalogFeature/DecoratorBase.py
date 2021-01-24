from balutils.CatalogFeature import CatalogFeature
from balutils import Catalog

class DecoratorBase(CatalogFeature):

    _catalogFeature: CatalogFeature = None

    def __init__(self, catalogFeature: CatalogFeature) -> None:
        self._catalogFeature = catalogFeature

    @property
    def component(self) -> None:
        return self._catalogFeature

    def applyTo(self, catalog: Catalog) -> None:
        return self._catalogFeature.operation(catalog)

