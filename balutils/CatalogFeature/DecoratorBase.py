from balutils.CatalogFeature import CatalogFeature

class DecoratorBase(CatalogFeature):

    _catalogFeature: CatalogFeature = None

    def __init__(self, catalogFeature: CatalogFeature) -> None:
        self._catalogFeature = catalogFeature

    @property
    def component(self) -> str:
        return self._catalogFeature

    def attach(self) -> None:
        return self._catalogFeature.operation()

