def dict_intersect(x,y): return {k:x[k] for k in set(x.keys()) & set(y.keys())}

def dict_union(x,y):
    d = {k:v for k,v in x.items()}
    d.update(y)
    return d