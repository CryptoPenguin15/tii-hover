# cython: language_level=3, boundscheck=False, wraparound=False
# distutils: libraries = m4ri
"""Bit-packed HOV degree-3/4 construction and in-place binary elimination."""
from sage.matrix.matrix_mod2_dense cimport Matrix_mod2_dense
from sage.libs.m4ri cimport mzd_read_bit, mzd_write_bit, mzd_t
from libc.stdlib cimport malloc, calloc, free
from libc.stdint cimport uint64_t
from cysignals.signals cimport sig_on, sig_off, sig_check
from sage.all import matrix, GF

cdef extern from "m4ri/m4ri.h":
    int mzd_echelonize(mzd_t *A, int full)

cdef int pair_index(int a, int b, int k):
    return a*(2*k-a-1)//2+b-a-1

def construct(Matrix_mod2_dense A, int p):
    cdef int k=A.nrows(), n=A.ncols(), a,b,c,d,j,t, stride
    cdef long col=0, N
    cdef unsigned char *v=<unsigned char*>malloc(k*n)
    cdef Matrix_mod2_dense E
    if v == NULL: raise MemoryError()
    try:
        for j in range(n):
            for a in range(k): v[j*k+a]=mzd_read_bit(A._entries,a,j)
        if p==3:
            N=k*(k-1)*(k-2)//6
            E=matrix(GF(2),n*k,N)
            for a in range(k):
                sig_check()
                for b in range(a+1,k):
                    for c in range(b+1,k):
                        for j in range(n):
                            t=j*k
                            if v[t+b] and v[t+c]: mzd_write_bit(E._entries,t+a,col,1)
                            if v[t+a] and v[t+c]: mzd_write_bit(E._entries,t+b,col,1)
                            if v[t+a] and v[t+b]: mzd_write_bit(E._entries,t+c,col,1)
                        col+=1
        elif p==4:
            N=(<long>k)*(k-1)*(k-2)*(k-3)//24
            stride=1+k*(k-1)//2
            E=matrix(GF(2),n*stride,N)
            for a in range(k):
                sig_check()
                for b in range(a+1,k):
                    for c in range(b+1,k):
                        for d in range(c+1,k):
                            for j in range(n):
                                t=j*k
                                if v[t+c] and v[t+d]: mzd_write_bit(E._entries,j*stride+1+pair_index(a,b,k),col,1)
                                if v[t+b] and v[t+d]: mzd_write_bit(E._entries,j*stride+1+pair_index(a,c,k),col,1)
                                if v[t+b] and v[t+c]: mzd_write_bit(E._entries,j*stride+1+pair_index(a,d,k),col,1)
                                if v[t+a] and v[t+d]: mzd_write_bit(E._entries,j*stride+1+pair_index(b,c,k),col,1)
                                if v[t+a] and v[t+c]: mzd_write_bit(E._entries,j*stride+1+pair_index(b,d,k),col,1)
                                if v[t+a] and v[t+b]: mzd_write_bit(E._entries,j*stride+1+pair_index(c,d,k),col,1)
                                if v[t+a] and v[t+b] and v[t+c] and v[t+d]: mzd_write_bit(E._entries,j*stride,col,1)
                            col+=1
        else: raise ValueError('Only p=3 and p=4 are implemented')
        return E
    finally: free(v)

def kernel_inplace(Matrix_mod2_dense E, int expected):
    """Consume E; avoid the extra full-matrix copy in a usual kernel call."""
    cdef int rank, r, p=0, f, j
    cdef int n=E.ncols()
    cdef Matrix_mod2_dense K
    sig_on()
    rank=mzd_echelonize(E._entries,1)
    sig_off()
    if n-rank != expected:
        raise ValueError('HOV kernel dimension %d, expected clean dimension %d'%(n-rank,expected))
    pivots=[]
    for r in range(rank):
        while p<n and not mzd_read_bit(E._entries,r,p): p+=1
        if p==n: raise ArithmeticError('Invalid echelon rank')
        pivots.append(p)
        p+=1
    pivot_set=set(pivots)
    frees=[j for j in range(n) if j not in pivot_set]
    K=matrix(GF(2),expected,n)
    for j,f in enumerate(frees):
        mzd_write_bit(K._entries,j,f,1)
        for r,p in enumerate(pivots):
            if mzd_read_bit(E._entries,r,f): mzd_write_bit(K._entries,j,p,1)
    return K

def replay(Matrix_mod2_dense A, Matrix_mod2_dense K, int p):
    """Independent streaming evaluation of all selected derivative constraints."""
    cdef int k=A.nrows(), n=A.ncols(), nb=K.nrows(), a,b,c,d,j,t,i,stride
    cdef long col=0, nr
    cdef uint64_t word
    cdef unsigned char *v=<unsigned char*>malloc(k*n)
    cdef uint64_t *out
    if nb>64: raise ValueError('Replay panel must have at most 64 rows')
    if p==3: stride=k
    elif p==4: stride=1+k*(k-1)//2
    else: raise ValueError('Only p=3 and p=4 are implemented')
    nr=n*stride
    out=<uint64_t*>calloc(nr,sizeof(uint64_t))
    if v==NULL or out==NULL:
        free(v); free(out); raise MemoryError()
    try:
        for j in range(n):
            for a in range(k): v[j*k+a]=mzd_read_bit(A._entries,a,j)
        for a in range(k):
            sig_check()
            for b in range(a+1,k):
                for c in range(b+1,k):
                    if p==3:
                        word=0
                        for i in range(nb): word|=(<uint64_t>mzd_read_bit(K._entries,i,col))<<i
                        if word:
                            for j in range(n):
                                t=j*k
                                if v[t+b] and v[t+c]: out[t+a]^=word
                                if v[t+a] and v[t+c]: out[t+b]^=word
                                if v[t+a] and v[t+b]: out[t+c]^=word
                        col+=1
                    else:
                        for d in range(c+1,k):
                            word=0
                            for i in range(nb): word|=(<uint64_t>mzd_read_bit(K._entries,i,col))<<i
                            if word:
                                for j in range(n):
                                    t=j*k
                                    if v[t+c] and v[t+d]: out[j*stride+1+pair_index(a,b,k)]^=word
                                    if v[t+b] and v[t+d]: out[j*stride+1+pair_index(a,c,k)]^=word
                                    if v[t+b] and v[t+c]: out[j*stride+1+pair_index(a,d,k)]^=word
                                    if v[t+a] and v[t+d]: out[j*stride+1+pair_index(b,c,k)]^=word
                                    if v[t+a] and v[t+c]: out[j*stride+1+pair_index(b,d,k)]^=word
                                    if v[t+a] and v[t+b]: out[j*stride+1+pair_index(c,d,k)]^=word
                                    if v[t+a] and v[t+b] and v[t+c] and v[t+d]: out[j*stride]^=word
                            col+=1
        return all(out[j]==0 for j in range(nr))
    finally:
        free(v); free(out)

def contraction(Matrix_mod2_dense K, int k, int p):
    """Raw tensor rows indexed by (p-1)-subsets; columns by (variable,basis)."""
    from itertools import combinations
    subsets=list(combinations(range(k),p-1))
    index={s:i for i,s in enumerate(subsets)}
    cdef int d=K.nrows(), j,r,i,row,col=0
    cdef Matrix_mod2_dense T=matrix(GF(2),len(subsets),k*d)
    for monomial in combinations(range(k),p):
        for r in monomial:
            row=index[tuple(x for x in monomial if x!=r)]
            for i in range(d):
                if mzd_read_bit(K._entries,i,col): mzd_write_bit(T._entries,row,r*d+i,1)
        col+=1
    return T

def store_equation_column(Matrix_mod2_dense E, int col, Matrix_mod2_dense P):
    cdef int a,b,n=P.ncols()
    for a in range(P.nrows()):
        for b in range(n):
            if mzd_read_bit(P._entries,a,b): mzd_write_bit(E._entries,a*n+b,col,1)
