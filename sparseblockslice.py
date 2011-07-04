class SparseBlocks(object):
    """ Slice object retunring chunks of input matrix defined by
    [chunksize] """
    def __init__(self, data, chunksize=[1]):
        self.data = data
        if sum(chunksize) != data.shape[0]:
            self.chunksize = [1]*data.shape[0]
            print 'Fallback default to chunks of size 1'
        else:
            self.chunksize = chunksize
    def __len__(self):
        length = len(self.chunksize)
        return length

    def _convert_slices(self, slices):
        newslices = []
        for axslice in slices:
            if isinstance(axslice, slice):
                start, stop = axslice.start, axslice.stop
                if axslice.start is not None:
                    start =  sum(self.chunksize[:axslice.start])
                if axslice.stop is not None:
                    stop = sum(self.chunksize[:axslice.stop+1])
                axslice = slice(start, stop, None)
            elif axslice is not None:
                axslice = slice(sum(self.chunksize[:axslice]), sum(self.chunksize[:axslice])+self.chunksize[axslice])
            newslices.append(axslice)
        return tuple(newslices)

    def __getitem__(self, item):
        item = self._convert_slices(item)
        return self.data.__getitem__(item)
    def __setitem__(self, item, value):
        item = self._convert_slices(item)
        return self.data.__setitem__(item, value)
