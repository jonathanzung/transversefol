using TransverseFol.Envelopes

swap = Envelopes.BasisChange([(0,0),(0,0)], [[1 0; 0 1], [0 1; -1 0]], [1,2])

# Upper Pareto front: (1/m, 3+3m) for m=1..20, (2,3), and their reflections in x==y
_nls_pts = vcat(
    [(Rational{Int}[1//m, 3+3m], nothing) for m in 1:20],
    [(Rational{Int}[2,   3    ], nothing)],
    [(Rational{Int}[3+3m, 1//m], nothing) for m in 1:20],
    [(Rational{Int}[3,   2    ], nothing)],
)
# Lower bound: 180-degree rotation of upper points, i.e. (a,b) -> (-a,-b)
_nls_pts_neg = [(Rational{Int}[-v[1], -v[2]], d) for (v, d) in _nls_pts]

El = Envelope{Lower}(_nls_pts_neg)
Eu = Envelope{Upper}(_nls_pts)

#multiview(["eLMkbcddddedde_2100_[(0,0),(0,0)]"], MM="m203", basis=swap, flows=[1,4], LS_envelope=[(El, Eu)], save_html=true)

multiview(["kvvLPQQkfghffijjijiaaaaaaabbbb_1020211100_[(0,0),(0,0),(3,-1)]","hLLLQkcegfeegghhhahabg_2010222_[(0,0),(0,0)]"], flows=[1,5], save_html=true)