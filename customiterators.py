from itertools import tee,izip_longest
def neighbour_zero(iterable):
    iterator = iter(iterable)
    prev = 0
    item = iterator.next()
    for next in iterator:
        yield (prev,item,next)
        prev = item
        item = next
    #yield (prev,item,0)

def pairwise(seq):
    a,b = tee(seq)
    b.next()
    return izip_longest(a,b)
