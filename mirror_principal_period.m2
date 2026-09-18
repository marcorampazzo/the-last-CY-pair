restart
needsPackage "NormalToricVarieties"
needsPackage "Polyhedra"


-- First, we define some custom methods
dot = (a,b) -> sum apply(#flatten entries a, i -> (flatten entries a)_i*(flatten entries b)_i)
encodeBinomial = f-> (transpose matrix exponents f)_0 - (transpose matrix exponents f)_1
createMonomial = (l, S) -> product apply(#l, i -> ((gens S)_i)^(l_i)) 
toricIdealFromTorusCharacters = (L, S)->(
    T := (coefficientRing S)[t_1..t_(rank target L), u_1..u_(rank target L)];
    T = T/ideal(toList(1..(rank target L))/(i -> t_i*u_i - 1)); -- coordinate ring of a torus
    parametrizingMonomials := apply(rank source L, j -> product apply(rank target L, i -> if L_(i,j) >= 0 then (t_(i+1))^(L_(i,j)) else (u_(i+1))^(-L_(i,j))));
    return ideal mingens ker map(T, S, parametrizingMonomials)
)


-- Next, we construct the ideals of the embeddings of the smooth Fano variety W-tilde and the toric variety W in P^11
K = ZZ/101
P11 = K[x_0..x_11]

M1 = random(P11^6,P11^{6:-1}); M1 = M1 - (transpose M1)
v1 = random(P11^6,P11^{1:-1})
idealW1 = ideal(M1*v1, pfaffians(6,M1))

M = matrix {
    {   0,  x_0,    0,    0,    0, -x_5},
    {-x_0,    0,  x_1,    0,    0,    0},
    {   0, -x_1,    0,  x_2,    0,    0},
    {   0,    0, -x_2,    0,  x_3,    0},
    {   0,    0,    0, -x_3,    0,  x_4},
    { x_5,    0,    0,    0, -x_4,    0}
}
v = transpose matrix{{x_6..x_11}}
idealW = ideal(M*v, pfaffians(6,M))

hilbertPolynomial idealW == hilbertPolynomial idealW1 -- true

integralRelns = sub(matrix transpose(((flatten entries mingens idealW)/(f -> encodeBinomial f))/(c -> flatten entries c)), QQ)
torusCharacters = transpose mingens ker transpose integralRelns 
toricIdealFromTorusCharacters(sub(torusCharacters, ZZ), P11) == idealW -- true


-- Next, we obtain the fan Σ that describes W as an abstract toric Fano variety
P = toSublattice convexHull torusCharacters -- Note that convexHull torusCharacters lives in an affine hyperplane in QQ^7. So, we apply "toSublattice" to it to see it in the smallest sub-lattice that contains it.
sigma = normalFan P
rank source rays sigma -- 18 = no. of rays of sigma

W = normalToricVariety sigma
isFano W -- true. So, W is Gorenstein Fano
facesAsCones(0, sigma)/(c -> rank source rays c == #(hilbertBasis c)) -- this is all true. So, W has at worst terminal singularities
rank picardGroup W -- 1


-- Next, we construct the three divisors H1, H2, H3 in the hyperplane class of W 
H1 = - toricDivisor(flatten entries ((fromCDivToWDiv W)*(transpose fromCDivToPic W)), W)
H2 = H1 + toricDivisor(flatten entries((rays sigma)^{0}), W)
H3 = H1 + toricDivisor(flatten entries((rays sigma)^{1}), W)
H1 + H2 + H3 == - toricDivisor W -- true

relnsAmongRays = mingens ker rays sigma
relnsAmongRaysAsCone = coneFromVData(relnsAmongRays|(-relnsAmongRays))
LSigma = intersection(relnsAmongRaysAsCone, posOrthant ambDim relnsAmongRaysAsCone)
gensLSigma = hilbertBasis LSigma


-- Next, we find the ideal of A(sigma) in C^18
R = K[t_0..t_(rank source rays sigma-1)]
conesSigma = join apply(1..dim sigma, d -> cones(d, sigma))
relns = flatten (conesSigma/(c->(
        temp := mingens ker((rays sigma)_c);
        apply(rank source temp, j -> product(apply(rank target temp, i -> if temp_(i,j) > 0 then (R_(c_i))^(temp_(i,j)) else 1)) - product(apply(rank target temp, i -> if temp_(i,j) < 0 then (R_(c_i))^(-temp_(i,j)) else 1)))
    )
))
idealASigma = saturate(ideal relns, product gens R)

generatingMonomials = gensLSigma/(g -> createMonomial(flatten entries g, R/idealASigma))
d = generatingMonomials / (m -> (degree m)_0); d = d/gcd(d)
gensLSigma/(f->dot(f,H1)) == d -- true


-- Next, we find the 12 intrinsic coordinates on L(sigma)
L = ZZ[l_0..l_17]
relnsAmongL = sub(reducedRowEchelonForm sub(hyperplanes LSigma, K), ZZ)
temp = mutableMatrix relnsAmongL;
toSub = {}
varsToEliminate = {}
for i in 0..((rank target relnsAmongL)-1) do (
    k = -1;
    for j in 0..((rank source relnsAmongL)-1) do if temp_(i,j) != 0 then (k = j; break);
    if k >= 0 then(
        varsToEliminate = append(varsToEliminate, L_k);
        temp_(i,k) = 0;
        toSub = append(toSub, L_k => - dot(temp^{i}, vars L));
    );
)
parameterizingVars = gens L;
varsToEliminate/(v -> parameterizingVars = delete(v, parameterizingVars))
inequalitiesDefiningLSigma = sub(halfspaces(LSigma)*(transpose vars L), toSub)
lDotH = sub(dot(vars L, H1), toSub) -- the intersection condition for nth coefficient is "lDotH" = n.