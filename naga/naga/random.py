from naga.naga.array import zeros


def random(shape, seed=1234):
    tmp = zeros(shape)
    tmp.data.init_prng()
    tmp.data.random(seed=seed)
    return tmp


def randn(shape, seed=1234):
    tmp = zeros(shape)
    tmp.data.init_prng()
    tmp.data.random(seed=seed, j='N')
    return tmp
