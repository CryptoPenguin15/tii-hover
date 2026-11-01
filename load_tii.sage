"""Load an original TII Track-2 public key into SageMath.

Usage: sage load_tii.sage data/pk_McEliece_100.txt
Importable Python implementation: from tii_cpu import read_public_key
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from tii_cpu import read_public_key

pk = read_public_key(sys.argv[1])
H, F = pk.sage_objects()
G = H.right_kernel().basis_matrix()
assert H * G.transpose() == 0
assert G.nrows() == pk.k
print('n={}, k={}, t={}, m={}'.format(pk.n, pk.k, pk.t, pk.m))
print('H:', H.nrows(), 'x', H.ncols())
print('G:', G.nrows(), 'x', G.ncols())
print('Field:', F)
print('SHA256:', pk.sha256)
