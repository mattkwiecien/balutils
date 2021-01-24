from balutils.CatalogFeature import SimpleCatalog
import numpy as np

class Detection(SimpleCatalog):
    """
    Adds detection catalog functionality to the catalog.
    """
    def applyTo(self, catalog: Catalog) -> None:
        
        self.parent.applyTo(catalog)

        '''
        Balrog stack versions 1.4 and below have a small bug that
        seems to duplicate exactly 1 object, so check for these
        '''
        unq, unq_idx, unq_cnt = np.unique(
            catalog._cat['bal_id'],
            return_inverse=True,
            return_counts=True
        )

        Nunq = len(unq)

        if Nunq != catalog.Nobjs:
            Ndups = catalog.Nobjs - Nunq
            dup_ids = unq[np.where(unq_cnt > 1)]
            print('Warning: Detection catalog has {} duplicate(s)!'.format(Ndups))
            print('Removing the following duplicates from detection catalog:')
            print(dup_ids)

            Nbefore = catalog.Nobjs
            for did in dup_ids:
                indx = np.where(catalog._cat['bal_id']==did)[0]

                L = len(indx)
                for i in range(L-1): # keep last one
                    catalog._cat.remove_row(indx[i])

            catalog.Nobjs = len(catalog._cat)
            assert catalog.Nobjs == (Nbefore - Ndups)

            print('{} duplicates removed, catalog size now {}'.format(Ndups, catalog.Nobjs))

        return

