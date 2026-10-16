import healpy as hp
import numpy as np

from fossil_sim.mapmaking import HitMaps


def test_hits_and_first_hit_time(tmp_path):
    maps = HitMaps((4, 8))
    vec = np.array(hp.pix2vec(8, np.array([0, 0, 5]))).T
    maps.add(vec, np.array([1, 2, 2], dtype=np.int8), np.array([30.0, 10.0, 20.0]))
    assert maps.total(8).sum() == 3
    assert maps.hits[8][1, 0] == 1 and maps.hits[8][2, 0] == 1
    assert maps.first_hit_time_s[8][0] == 10.0
    assert maps.coverage_fraction(8, np.array([15.0, 25.0]))[1] == 2 / hp.nside2npix(8)
    paths = maps.save(tmp_path)
    total = hp.read_map(paths[1], field=0, dtype=np.int64)
    assert total.sum() == 3
