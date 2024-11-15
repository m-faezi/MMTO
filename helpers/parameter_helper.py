import higra as hg
import numpy as np


def centroid(tree, size):
    emb = hg.EmbeddingGrid2d(size)
    coord = emb.lin2grid(np.arange(tree.num_leaves()))

    m = np.zeros((tree.num_leaves(), 3), dtype=np.float64)
    m[:, 0] = 1
    m[:, 1] = coord[:, 0]
    m[:, 2] = coord[:, 1]
    m = hg.accumulate_sequential(tree, m, hg.Accumulators.sum)
    m00 = m[:, 0]
    m10 = m[:, 1]
    m01 = m[:, 2]

    xmean = m10 / m00
    ymean = m01 / m00

    return xmean[tree.num_leaves():], ymean[tree.num_leaves():]


def weighted_centroid(tree, size, image):
    image -= max(np.min(image), 0)
    im = np.exp(image)
    emb = hg.EmbeddingGrid2d(size)
    coord = emb.lin2grid(np.arange(tree.num_leaves()))

    m = np.zeros((tree.num_leaves(), 3), dtype=np.float64)
    m[:, 0] = 1
    m[:, 1] = coord[:, 0]
    m[:, 2] = coord[:, 1]
    m = np.array([[
        im[int(i[1]), int(i[2])],
        i[1] * im[int(i[1]), int(i[2])],
        i[2] * im[int(i[1]), int(i[2])]
    ] for i in m])
    m = hg.accumulate_sequential(tree, m, hg.Accumulators.sum)
    m00 = m[:, 0]
    m10 = m[:, 1]
    m01 = m[:, 2]

    xmean = m10 / m00
    ymean = m01 / m00

    return xmean[tree.num_leaves():], ymean[tree.num_leaves():]


def total_flux(tree, size, image):
    emb = hg.EmbeddingGrid2d(size)
    coord = emb.lin2grid(np.arange(tree.num_leaves()))

    m = np.zeros((tree.num_leaves(), 3), dtype=np.float64)
    m[:, 0] = 1
    m[:, 1] = coord[:, 0]
    m[:, 2] = coord[:, 1]
    m = np.array([[
        image[int(i[1]), int(i[2])],
        i[1] * image[int(i[1]), int(i[2])],
        i[2] * image[int(i[1]), int(i[2])]
    ] for i in m])
    m = hg.accumulate_sequential(tree, m, hg.Accumulators.sum)
    t_f = m[:, 0]

    return t_f[tree.num_leaves():]

