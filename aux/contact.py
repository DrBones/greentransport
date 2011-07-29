import functools
import scipy as sc
class Contact(sc.ndarray):
    """Simple wrapper class to add SO attribute, otherwise
    behaves like ndarray, all below is probably overkill, less
    would suffice but i dont't know enough about all this"""
    def __new__(cls, input_array, SO=None):
        obj =  sc.asarray(input_array).view(cls)
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.SO = getattr(obj, 'SO', None)
