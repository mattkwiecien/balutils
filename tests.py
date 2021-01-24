from balutils.CatalogFactory import CatalogFactory
from balutils import Catalog

def test_create_fits():
    cat = CatalogFactory.Create('testfile.nm', ['fits'])

def main():
    test_create_fits()
    # test_create_detection()
    # test_create_gold_default()
    # test_create_gold_sof()
    # test_create_gold_mof()

if __name__ == '__main__':
    main()
