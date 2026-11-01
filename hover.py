#!/usr/bin/env python3
"""CPU HOV/HOVER recovery following Saarinen, ePrint 2026/1778.

Binary HOV degrees 3 and 4, certified field selector, degree-3 scroll
completion, and full block-omitting-window completion. All keys must pass the
independent binary Goppa verifier before output. This is a new implementation;
timings from Saarinen's implementation are not performance claims for this one.
"""
import argparse
import gc
import itertools
import json
import math
from pathlib import Path
import random
import sys
import time

from sage.all import GF, PolynomialRing, matrix, vector, identity_matrix
from sage.misc.cython import cython_import

from tii_cpu import (PublicKey, read_public_key, binary_rref, parity_rows,
                     verify_key, BinaryField)
from recover_key import integer_element, int_coefficients, write_key, atomic_json

NATIVE = None


def native():
    global NATIVE
    if NATIVE is None:
        NATIVE = cython_import(str(Path(__file__).with_name('hov_native.pyx')),verbose=0,use_cache=True)
    return NATIVE


def say(message):
    print(time.strftime('%H:%M:%S')+' '+str(message),flush=True)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def save_kernel(path,pk,pr,K):
    rows=[]
    for row in K:
        packed=bytearray((K.ncols()+7)//8)
        for i,v in enumerate(row):
            if v: packed[i//8]|=1<<(i%8)
        rows.append(packed.hex())
    atomic_json(path,dict(public_key_sha256=pk.sha256,profile=pr,
                          ncols=K.ncols(),rows_lsb_hex=rows))


def read_kernel(path,pk,pr):
    obj=json.loads(path.read_text())
    require(obj['public_key_sha256']==pk.sha256 and obj['profile']==pr,
            'Checkpoint provenance mismatch')
    require(obj['ncols']==pr['columns'] and len(obj['rows_lsb_hex'])==pk.m,
            'Checkpoint dimensions are invalid')
    rows=[]
    for text in obj['rows_lsb_hex']:
        raw=bytes.fromhex(text)
        require(len(raw)==(pr['columns']+7)//8,'Invalid packed checkpoint row')
        rows.append([(raw[i//8]>>(i%8))&1 for i in range(pr['columns'])])
    return matrix(GF(2),rows)


def profile(pk,p):
    # Integer arithmetic avoids a rounding error at powers of two in eq. (3).
    e=1
    while (p-1)*(2**e) <= pk.t-1:
        e+=1
    # e = 1 + floor(log2((t-1)/(p-1))).
    s=(e+2)*pk.t-(p-1)*2**e-p
    kappa=len(pk.rows)-s
    require(e>=1 and 0<=s<len(pk.rows),'Invalid HOV parameter profile')
    a=pk.m-e-2
    require(a>0,'Empty HOV window')
    columns=math.comb(kappa,p)
    rows=pk.k*sum(math.comb(kappa,p-2**ell)
                 for ell in range(1,p.bit_length()) if 2**ell<=p)
    return dict(p=p,e=e,s=s,kappa=kappa,window=a,
                rank_bound=(pk.m-1)*pk.t-s,rows=rows,columns=columns,
                packed_matrix_GiB=rows*((columns+63)//64)*8/2**30)


def choose_degree(pk,requested):
    if requested!='auto':
        return int(requested)
    for p in (3,4):
        pr=profile(pk,p)
        if pr['rows']>=pr['columns']-pk.m:
            return p
    raise ValueError('Neither degree 3 nor degree 4 has a reachable clean-kernel profile')


def shortening(H,pr,seed):
    r,n=H.dimensions()
    order=list(range(n))
    random.Random(seed).shuffle(order)
    reduced=H.matrix_from_columns(order).echelon_form()
    pivots=list(reduced.pivots())
    require(len(pivots)==r,'Public parity-check rank changed')
    nonpivots=[i for i in range(n) if i not in set(pivots)]
    perm=pivots+nonpivots
    systematic=reduced.matrix_from_columns(perm)
    require(systematic[:,:r]==identity_matrix(GF(2),r),'Failed systematic form')
    original_indices=[order[i] for i in perm]
    short=systematic[pr['s']:,:]
    lifted=matrix(GF(2),pr['kappa'],n)
    for j,original in enumerate(original_indices):
        lifted.set_column(original,short.column(j))
    A=short[:,r:]
    require(A.ncols()==n-r,'Wrong nonpivot-column count')
    return A,lifted


def hov_kernel(A,pr,m):
    say('Constructing HOV matrix {} x {} ({:.3f} GiB packed)'.format(
        pr['rows'],pr['columns'],pr['packed_matrix_GiB']))
    E=native().construct(A,pr['p'])
    say('Computing HOV kernel with in-place M4RI elimination')
    K=native().kernel_inplace(E,m)
    del E
    gc.collect()
    require(native().replay(A,K,pr['p']),'HOV kernel replay failed')
    say('Full selected-derivative replay passed; kernel dimension {}'.format(K.nrows()))
    return K


def random_control(pr,nonpivots,seed):
    """Matched random systematic dual; require zero HOV kernel, as in the paper."""
    from sage.all import set_random_seed, random_matrix
    set_random_seed(seed)
    A=random_matrix(GF(2),pr['kappa'],nonpivots)
    say('Checking matched random-code profile has zero HOV kernel')
    E=native().construct(A,pr['p'])
    K=native().kernel_inplace(E,0)
    del E
    gc.collect()
    require(K.nrows()==0,'Matched random-code kernel is nonzero')
    say('Matched random-code control passed')


def catalecticant(coefficients,k,p,F):
    subsets=list(itertools.combinations(range(k),p-1))
    index={s:i for i,s in enumerate(subsets)}
    C=matrix(F,k,len(subsets))
    for monomial,c in zip(itertools.combinations(range(k),p),coefficients):
        if c:
            for r in monomial:
                C[r,index[tuple(j for j in monomial if j!=r)]]=c
    return C


def selector(A,K,lifted,pk,pr):
    require(K.nrows()==pk.m,'HOVER requires a clean dimension-m kernel')
    raw=native().contraction(K,pr['kappa'],pr['p'])
    T=raw.row_space().basis_matrix()
    del raw
    Z=T.right_kernel_matrix()
    d=pk.m
    say('Contraction rank {}; compatible-algebra system {} x {}'.format(
        T.nrows(),T.nrows()*Z.nrows(),d*d))
    equations=matrix(GF(2),T.nrows()*Z.nrows(),d*d)
    blocksT=[T.matrix_from_columns(list(range(i,T.ncols(),d))) for i in range(d)]
    blocksZ=[Z.matrix_from_columns(list(range(i,Z.ncols(),d))).transpose() for i in range(d)]
    for i in range(d):
        for j in range(d):
            native().store_equation_column(equations,i*d+j,blocksT[i]*blocksZ[j])
    algebra=equations.right_kernel_matrix()
    require(algebra.nrows()==d,'Compatible algebra dimension {}, expected {}'.format(algebra.nrows(),d))
    I=identity_matrix(GF(2),d)
    space=algebra.row_space()
    require(vector(GF(2),I.list()) in space,'Identity not in compatible algebra')
    Bs=[matrix(GF(2),d,d,list(row)) for row in algebra]
    for B in Bs:
        for C in Bs:
            require(vector(GF(2),(B*C).list()) in space,'Compatible algebra is not closed')
    generator=None
    rng=random.Random(0)
    for trial in range(256):
        B=Bs[trial] if trial<d else sum((b for b in Bs if rng.randrange(2)),matrix(GF(2),d,d))
        if B.charpoly().is_irreducible():
            generator=B
            break
    require(generator is not None,'No irreducible field generator found')
    _,F=pk.sage_objects()
    roots=generator.charpoly().change_ring(F).roots(multiplicities=False)
    require(len(roots)==d,'Field generator failed to split in the TII field')
    root=min(roots,key=integer_element)
    eigenvectors=[]
    B=generator.change_ring(F)
    for i in range(d):
        eig=(B-root**(2**i)*identity_matrix(F,d)).right_kernel_matrix()
        require(eig.nrows()==1,'An eigenspace is not one dimensional')
        v=eig.row(0)
        first=next(c for c in v if c)
        eigenvectors.append(v/first)
    for i in range(d):
        squared=vector(F,[x*x for x in eigenvectors[i]])
        require(squared==eigenvectors[(i+1)%d],'Eigenlines do not form one Frobenius cycle')
    say('Field certificate passed: algebra GF(2^{})'.format(d))
    KF=K.change_ring(F)
    cats=[]
    directions=[]
    cat_rank_sum=0
    for i,v in enumerate(eigenvectors):
        C=catalecticant(v*KF,pr['kappa'],pr['p'],F)
        U=C.column_space().basis_matrix()
        cat_rank_sum+=U.nrows()
        directions.append(U*lifted.change_ring(F))
        cats.append(C)
        say('Direction {}: catalecticant rank {}'.format(i,U.nrows()))
    require(cat_rank_sum==T.nrows(),'Eigenline images fail the direct-sum rank certificate')
    kernels=[C.left_kernel() for C in cats]
    require(all(kernels[i]!=kernels[j] for i in range(d) for j in range(i)),
            'Derivative kernels do not separate all directions')
    windows=[]
    for i in range(d):
        C=sum((cats[(i+j)%d] for j in range(pr['window'])),matrix(F,cats[0].nrows(),cats[0].ncols()))
        U=C.column_space().basis_matrix()
        require(U.nrows()<=pr['rank_bound'],'Window fails the block-omission rank bound')
        windows.append(U*lifted.change_ring(F))
    say('All Frobenius-window rank gates passed')
    return dict(directions=directions,windows=windows,profile=pr,
                contraction_rank=T.nrows(),algebra_dimension=d,
                direction_ranks=[C.rank() for C in cats],
                window_ranks=[W.rank() for W in windows],
                field_generator=[int(x) for x in generator.list()],
                irreducible_characteristic_polynomial=[int(x) for x in generator.charpoly().list()])


def restrict_quadrics(Q,pairs,U):
    """Coefficient matrix after substituting x=z*U into homogeneous quadrics."""
    F=Q.base_ring()
    r=U.nrows()
    small=list(itertools.combinations_with_replacement(range(r),2))
    transform=matrix(F,len(pairs),len(small))
    for q,(i,j) in enumerate(pairs):
        for s,(a,b) in enumerate(small):
            transform[q,s]=(U[a,i]*U[b,j] if a==b else
                            U[a,i]*U[b,j]+U[b,i]*U[a,j])
    return Q*transform,small


def ruling_at_point(Q,pairs,point,p=3):
    """Recover the scroll's P^2 ruling from its tangent quadrics in char 2.

    On the (p+1)-dimensional affine tangent space the restricted quadratic
    ideal has a unique square direction l^2. Its square root l cuts out the
    p-dimensional ruling. Every resulting plane is checked against all
    public quadrics before it is used.
    """
    F=Q.base_ring()
    J=matrix(F,Q.nrows(),p*p)
    for s,(i,j) in enumerate(pairs):
        if i!=j:
            for q in range(Q.nrows()):
                J[q,i]+=Q[q,s]*point[j]
                J[q,j]+=Q[q,s]*point[i]
    U=J.right_kernel_matrix()
    require(U.nrows()==p+1,'Scroll tangent dimension does not equal p+1')
    restricted,small=restrict_quadrics(Q,pairs,U)
    off=[i for i,(a,b) in enumerate(small) if a!=b]
    diag=[i for i,(a,b) in enumerate(small) if a==b]
    combinations=restricted.matrix_from_columns(off).left_kernel_matrix()
    squares=(combinations*restricted).matrix_from_columns(diag).row_space().basis_matrix()
    require(squares.nrows()==1,'Scroll tangent square direction is not unique')
    ell=matrix(F,1,p+1,[v.sqrt() for v in squares.row(0)])
    plane=ell.right_kernel_matrix()*U
    check,_=restrict_quadrics(Q,pairs,plane)
    require(check.is_zero(),'Recovered ruling does not lie in the quadratic ideal')
    require(plane.rank()==p and vector(F,point) in plane.row_space(),'Ruling does not contain its point')
    return plane


def scroll_support(V,p=3):
    F=V.base_ring()
    visible=[j for j in range(V.ncols()) if not V.column(j).is_zero()]
    ambient=p*p
    require(V.nrows()==ambient and V.rank()==ambient,
            'Direction rank does not equal p^2 for the scroll completion')
    pairs=list(itertools.combinations_with_replacement(range(ambient),2))
    hilbert=(2*p-1)*p*(p+1)//2
    quadrics=len(pairs)-hilbert
    require(len(visible)>=hilbert,'Too few visible points for scroll recognition')
    evals=matrix(F,[[V[i,j]*V[k,j] for i,k in pairs] for j in visible])
    Q=evals.right_kernel_matrix()
    require(Q.nrows()==quadrics,'Quadratic ideal dimension does not match the scroll')
    say('Scroll quadratic ideal recognized; recovering labeled rulings')
    triples=list(itertools.combinations(range(ambient),p))
    curve_degree=p*(p-1)
    planes=[]
    for j in visible[:curve_degree+5]:
        planes.append(ruling_at_point(Q,pairs,V.column(j),p))
    samples=matrix(F,[[P.matrix_from_columns(t).det() for t in triples] for P in planes])
    require(samples.rank()==curve_degree+1,'Ruling Pluecker rank does not match a rational-normal curve')
    keep=list(samples.pivots())
    coords=[triples[i] for i in keep]
    columns=[]
    for index,j in enumerate(visible):
        P=planes[index] if index<len(planes) else ruling_at_point(Q,pairs,V.column(j),p)
        columns.append([P.matrix_from_columns(t).det() for t in coords])
        if index and index%200==0:
            say('Recovered {} / {} ruling labels'.format(index,len(visible)))
    grs=matrix(F,columns).transpose()
    support=grs_support(grs)
    return visible,support


def grs_support(C):
    """Sidelnikov--Shestakov cross ratios, followed by exact GRS verification."""
    F=C.base_ring()
    reduced=C.echelon_form()
    k=reduced.rank()
    reduced=reduced[:k,:]
    require(k>=3 and reduced.ncols()>=k+2,'GRS extraction has insufficient dimension or length')
    piv=list(reduced.pivots())
    non=[j for j in range(C.ncols()) if j not in set(piv)]
    a,b=non[:2]
    require(all(reduced[i,j] for i in range(k) for j in non),'GRS systematic minor contains zero')
    scale=reduced[0,a]/reduced[1,a]
    beta=[None]*C.ncols()
    beta[piv[0]]=None  # projective infinity
    beta[piv[1]]=F.zero()
    for j in non:
        beta[j]=reduced[0,j]/reduced[1,j]/scale
    for i in range(2,k):
        ratio=(reduced[i,a]/reduced[1,a])/(reduced[i,b]/reduced[1,b])
        den=beta[a]+ratio*beta[b]
        require(bool(den),'GRS pivot support extraction is singular')
        beta[piv[i]]=beta[a]*beta[b]*(1+ratio)/den
    finite={x for x in beta if x is not None}
    require(len(finite)==C.ncols()-1,'Recovered projective GRS support is not distinct')
    pole=next((F.from_integer(i) for i in range(F.cardinality()) if F.from_integer(i) not in finite),None)
    require(pole is not None,'No affine chart avoids the recovered support')
    alpha=[F.zero() if x is None else 1/(x+pole) for x in beta]
    W=matrix(F,[[x**i for x in alpha] for i in range(k)])
    W=W.matrix_from_columns(piv).inverse()*W
    multipliers=[F.zero()]*C.ncols()
    multipliers[piv[0]]=F.one()
    for j in non:
        multipliers[j]=reduced[0,j]/W[0,j]
    for i in range(1,k):
        multipliers[piv[i]]=W[i,a]*multipliers[a]/reduced[i,a]
    require(all(reduced[i,j]*multipliers[piv[i]]==W[i,j]*multipliers[j]
                for i in range(k) for j in range(C.ncols())),
            'Recovered support and multipliers fail exact GRS reconstruction')
    return alpha


def finish_partial(H,pk,visible,beta):
    """Recover the Goppa polynomial and all omitted labels from partial support."""
    F=beta[0].parent()
    R=PolynomialRing(F,'X')
    X=R.gen()
    restricted=H.matrix_from_columns(visible)
    require(restricted.rank()==len(pk.rows),'Visible public matrix is rank deficient')
    even=restricted.stack(matrix(GF(2),1,len(visible),[1]*len(visible))).right_kernel_matrix()
    require(even.nrows()>0,'No even-weight shortened codewords')
    common=R.zero()
    for row in even:
        product=R.one()
        for j,c in enumerate(row):
            if c: product*=X+beta[j]
        derivative=product.derivative()
        if derivative:
            common=common.gcd(derivative)
            require(common.degree()>=2*pk.t,'Even-codeword gcd has lost the Goppa factor')
            if common.degree()==2*pk.t: break
    require(bool(common),'Even-codeword derivatives all vanish')
    root=R([common[2*i].sqrt() for i in range(common.degree()//2+1)])
    factors=[f.monic() for f,_ in root.factor() if f.degree()==pk.t]
    require(len(factors)==1,'Goppa factor is not uniquely determined')
    g=factors[0]
    odd=next((row for row in restricted.right_kernel_matrix() if sum(map(int,row))%2),None)
    pole=None  # None denotes infinity.
    if odd is not None:
        product=R.one()
        for j,c in enumerate(odd):
            if c: product*=X+beta[j]
        derivative=product.derivative()%g
        if derivative:
            residue=(X+product*derivative.inverse_mod(g))%g
            require(residue.degree()<=0,'Projective pole is not a field constant')
            pole=residue[0]
    if pole is not None:
        require(pole not in beta,'Projective pole lies on the visible support')
        alpha=[1/(x+pole) for x in beta]
        g=sum((g[i]*X**(pk.t-i)*(1+pole*X)**i for i in range(pk.t+1)),R.zero()).monic()
    else:
        alpha=beta
    require(g.is_irreducible(),'Recovered Goppa polynomial is not irreducible')
    field=pk.field
    gg=int_coefficients(g)
    known=[integer_element(x) for x in alpha]
    recovered=matrix(GF(2),len(pk.rows),len(known),
                     [[(row>>j)&1 for j in range(len(known))]
                      for row in parity_rows(known,gg,field)])
    piv=list(restricted.pivots())
    transform=recovered.matrix_from_columns(piv)*restricted.matrix_from_columns(piv).inverse()
    require(transform*restricted==recovered,'Partial support fails public row-operation consistency')
    full=transform*H
    support=[]
    for j in range(pk.n):
        y=sum(int(full[bit,j])<<bit for bit in range(pk.m))
        xy=sum(int(full[pk.m+bit,j])<<bit for bit in range(pk.m))
        require(y!=0,'A completed multiplier is zero')
        support.append(field.mul(xy,field.inverse(y)))
    verify_key(pk,gg,support)
    return gg,support


def frobenius_space(M):
    return matrix(M.base_ring(),M.nrows(),M.ncols(),[x*x for x in M.list()]).row_space()


def window_completion(H,pk,records):
    F=records[0]['windows'][0].base_ring()
    target=(pk.m-1)*pk.t
    # Global Frobenius phase is irrelevant. Fix the first shortening's phase.
    states=[records[0]['windows'][0].row_space()]
    for rec in records[1:]:
        updated=[]
        for accumulated in states:
            for W in rec['windows']:
                combined=accumulated+W.row_space()
                if combined.dimension()<=target and combined not in updated:
                    updated.append(combined)
        states=updated
        say('Window phase states surviving exact dimension gate: {}'.format(len(states)))
        for W in states:
            if W.dimension()!=target:
                continue
            extracted=W
            shifted=W
            valid=True
            for step in range(1,pk.m-1):
                shifted=frobenius_space(shifted.basis_matrix())
                extracted=extracted.intersection(shifted)
                if extracted.dimension()!=(pk.m-1-step)*pk.t:
                    valid=False
                    break
            if valid:
                try:
                    alpha=grs_support(extracted.basis_matrix())
                    return finish_partial(H,pk,list(range(pk.n)),alpha)
                except ValueError as error:
                    say('Window completion candidate refused: {}'.format(error))
    return None


def recover(pk,p,trials,seed,checkpoint,max_gib,check_random=True):
    H,F=pk.sage_objects()
    pr=profile(pk,p)
    say(json.dumps(pr,sort_keys=True))
    require(pr['rows']>=pr['columns']-pk.m,'Profile has too many unconstrained monomials for a clean kernel')
    require(pr['packed_matrix_GiB']<=max_gib,'HOV matrix exceeds --max-matrix-gib')
    if check_random:
        random_control(pr,pk.k,seed+1234567)
    records=[]
    for trial in range(trials):
        say('Public shortening trial {}'.format(trial))
        A,lifted=shortening(H,pr,seed+trial)
        cache=checkpoint/('hov_{}_p{}_seed{}.json'.format(pk.sha256[:16],p,seed+trial)) if checkpoint else None
        if cache is not None and cache.exists():
            K=read_kernel(cache,pk,pr)
            require(native().replay(A,K,p),'Loaded HOV kernel fails full public replay')
            say('Loaded and replayed public-derived kernel checkpoint')
        else:
            try:
                K=hov_kernel(A,pr,pk.m)
            except ValueError as error:
                say('Shortening refused: {}'.format(error))
                continue
            if cache is not None:
                checkpoint.mkdir(parents=True,exist_ok=True)
                save_kernel(cache,pk,pr,K)
        try:
            rec=selector(A,K,lifted,pk,pr)
        except ValueError as error:
            say('Selector refused: {}'.format(error))
            continue
        records.append(rec)
        if cache is not None:
            certificate={name:value for name,value in rec.items() if name not in ('directions','windows')}
            certificate.update(public_key_sha256=pk.sha256,seed=seed+trial,
                               random_control_checked=check_random)
            atomic_json(cache.with_suffix('.certificate.json'),certificate)
        for index,V in enumerate(rec['directions']):
            try:
                visible=[j for j in range(V.ncols()) if not V.column(j).is_zero()]
                try:
                    alpha=grs_support(V.matrix_from_columns(visible))
                    say('Direction is already an exactly reconstructed GRS code')
                except ValueError:
                    visible,alpha=scroll_support(V,p)
                key=finish_partial(H,pk,visible,alpha)
                say('One-direction completion passed full public-key verification')
                return key
            except ValueError as error:
                say('Direction {} completion refused: {}'.format(index,error))
        if len(records)>1:
            key=window_completion(H,pk,records)
            if key is not None:
                return key
    raise ValueError('No certified completion after the requested public shortenings')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('public_key',type=Path)
    parser.add_argument('--degree',choices=['auto','3','4'],default='auto')
    parser.add_argument('--trials',type=int,default=3)
    parser.add_argument('--seed',type=int,default=0)
    parser.add_argument('--output',type=Path,default=Path('hover_key.json'))
    parser.add_argument('--checkpoint-dir',type=Path)
    parser.add_argument('--max-matrix-gib',type=float,default=6.0)
    parser.add_argument('--profile-only',action='store_true')
    parser.add_argument('--skip-random-control',action='store_true',
                        help='Skip the matched random-code check (a diagnostic option)')
    args=parser.parse_args()
    pk=read_public_key(args.public_key)
    p=choose_degree(pk,args.degree)
    if args.profile_only:
        print(json.dumps(profile(pk,p),indent=2))
        return 0
    try:
        key=recover(pk,p,args.trials,args.seed,args.checkpoint_dir,args.max_matrix_gib,
                    check_random=not args.skip_random_control)
        write_key(args.output,pk,key,'HOV/HOVER degree {}'.format(p))
    except ValueError as error:
        say('INCONCLUSIVE: '+str(error))
        return 2
    return 0


if __name__=='__main__':
    sys.exit(main())
