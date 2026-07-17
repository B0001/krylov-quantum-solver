// ============================================================
// KNOWLEDGE GRAPH v2 — Computational Complexity of the
// Electronic Structure Problem
// arXiv-verified edition, frontier-checked 2026-07-16
// ============================================================
//
// CHANGELOG v1 -> v2
//   [FIXED] es_fixedbasis: HARD_FOR -> COMPLETE_FOR (QMA).
//     arXiv:2103.08215 proves full QMA-COMPLETENESS. v1's note
//     ("Coulomb tensor fixed, hardness via one-body terms") was
//     WRONG: the construction sets V = 0 and encodes the hard
//     instance entirely in the CHOICE OF GAUSSIAN ORBITALS.
//     Verified against the paper's Definition 1 and Theorem 1.
//   [ADDED] Full reduction chain: QMA circuit -> LH5 -> ... ->
//     antiferromagnetic Heisenberg (Piddock-Montanaro) ->
//     Fermi-Hubbard (O'Gorman et al. Thm 3) -> ESFBS (Thm 1).
//   [FIXED] LESD/Hartree-Fock: NP-completeness for electronic-
//     structure Hamiltonians specifically is O'Gorman et al.
//     Thm 2; Schuch-Verstraete showed NP-hardness for generic
//     quartic fermionic Hamiltonians.
//   [ADDED] mol_gsee is no longer edge-isolated: its OPEN
//     status is now CITABLE — arXiv:2103.08215 explicitly poses
//     the nuclei-only-potential variant as an open problem.
//   [ADDED] Guided-LH refinements 2022-2025, incl. arXiv:
//     2207.10097, 2207.10250, 2302.11578, 2411.16163,
//     2509.25815, 2509.25829.
//   [ADDED] New axis: conditional SAMPLING-hardness for
//     chemistry circuits (2503.21041, 2504.12893) — distinct
//     from decision-problem hardness.
//   [ADDED] Evidence-debate frontier: Lee et al. 2208.02199;
//     FeMoco solved classically to chemical accuracy
//     (2601.04621, Jan 2026); Fe4S4 Blackwell-DMRG benchmark
//     (2603.28648, Mar 2026).
//   [VERDICT] Frontier sweep 2026-07-16: no paper found that
//     resolves the exact complexity of Mol-GSEE. Open core
//     remains open.
//
// CITATION PROVENANCE
//   Edges marked verified:"2026-07-16" were checked this
//   session against the arXiv full text or alphaXiv discovery.
//   Other arXiv IDs are standard bibliography identifiers for
//   classic results (pre-2020), included for convenience but
//   not re-fetched this session.
//
// USAGE
//   One Cypher query; all CREATE clauses share variable scope.
//   Paste whole into Neo4j Browser, or session.run(file). For
//   statement-by-statement execution, convert relationship
//   clauses to MATCH-by-id form; every node carries unique `id`.
//
// STATUS TAG LEGEND
//   THM proven theorem | DEF by definition | FACT standard
//   result | COROLLARY immediate consequence | EMPIRICAL
//   observed, unproven | CONJ conjectured | OPEN unresolved |
//   DESIGN engineering rationale | META about this graph
//
// MODELING NOTES
//   - Complexity attaches to PROBLEMS; Problem is first-class.
//   - Definitional content lives in node properties, not edges.
//   - "Amorphous boundary conditions": no established referent;
//     kept as flagged MetaNote only.
//   - This graph models an argument's structure; it is not
//     itself a proof.
// ============================================================


// ---------- A. OPERATORS, SPACES, TRANSFORMS ----------
CREATE (h_el:Operator {id:"h_el", name:"Electronic Hamiltonian",
  form:"Sum h_pq a†_p a_q + (1/2) Sum g_pqrs a†_p a†_q a_r a_s (two-body fermionic)"})
CREATE (fock:Space {id:"fock", name:"Fock Space (finite orbital basis)",
  dimension:"C(M,N_alpha)*C(M,N_beta)", scaling:"exponential in N, for M = Theta(N)"})
CREATE (config_l2:Space {id:"config_l2", name:"Antisymmetric L2(R^3N) continuum space"})
CREATE (qubit_h:Operator {id:"qubit_h", name:"Qubit-encoded Hamiltonian"})
CREATE (basis_projection:Process {id:"basis_projection", name:"Finite-basis projection"})
CREATE (basis_incompleteness:Error {id:"basis_incompleteness", name:"Basis-set incompleteness error"})

CREATE (h_el)-[:ACTS_ON {status:"DEF"}]->(fock)
CREATE (fock)-[:PROJECTS {status:"DEF",
  note:"finite-basis Fock space is a projection of the continuum problem"}]->(config_l2)
CREATE (basis_projection)-[:INTRODUCES {status:"FACT"}]->(basis_incompleteness)
CREATE (h_el)-[:MAPS_POLY_TO {status:"THM",
  note:"Jordan-Wigner / Bravyi-Kitaev encoding"}]->(qubit_h)


// ---------- B. PROMISE PROBLEMS + COMPLEXITY CLASSES ----------
CREATE (lh5:Problem {id:"lh5", name:"k-Local Hamiltonian, k=5"})
CREATE (lh2:Problem {id:"lh2", name:"k-Local Hamiltonian, k=2"})
CREATE (lh2_2d:Problem {id:"lh2_2d", name:"2-Local Hamiltonian, 2D qubit lattice"})
CREATE (lh_1d_d12:Problem {id:"lh_1d_d12", name:"Local Hamiltonian, 1D nearest-neighbor, local dim 12"})
CREATE (lh_1d_ti:Problem {id:"lh_1d_ti", name:"Local Hamiltonian, 1D translationally-invariant"})
CREATE (lh_stoq:Problem {id:"lh_stoq", name:"Local Hamiltonian, stoquastic"})
CREATE (lh_stoq_ff:Problem {id:"lh_stoq_ff", name:"Local Hamiltonian, stoquastic + frustration-free"})
CREATE (lh_constgap:Problem {id:"lh_constgap", name:"Local Hamiltonian, constant relative promise gap",
  definition:"the quantum-PCP regime"})
CREATE (heisenberg_afm:Problem {id:"heisenberg_afm",
  name:"Antiferromagnetic Heisenberg Hamiltonian",
  definition:"ground-state energy of Sum kappa_ij (XX+YY+ZZ)_ij, kappa_ij >= 0 poly-bounded"})
CREATE (fermi_hubbard:Problem {id:"fermi_hubbard",
  name:"Fermi-Hubbard, uniform onsite repulsion",
  definition:"ground energy in n-particle subspace; uniform u0, arbitrary weighted graph, |t_ij| poly-bounded, all hoppings may share one sign; no magnetic field"})
CREATE (guidedlh_poly:Problem {id:"guidedlh_poly", name:"Guided Local Hamiltonian, 1/poly precision",
  definition:"LH plus classical guiding state u with |<u|psi_0>| >= 1/poly(n); output precision 1/poly(n)"})
CREATE (guidedlh_const:Problem {id:"guidedlh_const", name:"Guided Local Hamiltonian, constant precision",
  definition:"as guidedlh_poly but constant output precision"})
CREATE (guidedlh_stoq:Problem {id:"guidedlh_stoq", name:"Guided Local Hamiltonian, stoquastic"})
CREATE (guidable_lh:Problem {id:"guidable_lh", name:"Guidable Local Hamiltonian",
  definition:"Merlinized variant: existence of a good guiding state is promised, not given"})
CREATE (es_v:Problem {id:"es_v", name:"Electronic Structure, tunable external potential",
  definition:"ground-state energy of electrons in tunable external (magnetic-field-augmented) potential, promise gap >= 1/poly"})
CREATE (es_fixedbasis:Problem {id:"es_fixedbasis",
  name:"Electronic Structure in a Fixed Basis Set (ESFBS)",
  definition:"instance = (V, eta electrons, basis phi_1..phi_n, thresholds a<b, b-a >= 1/poly(eta)); integrals t_ij, v_ij, u_ijkl efficiently computable from concise basis description; decide E0 <= a vs >= b in the eta-electron span",
  key_fact:"hard instances use V = 0 with the QMA-hard content encoded entirely in the CHOICE of Gaussian-superposition orbitals; arbitrary total spin is essential to the construction",
  verified:"2026-07-16"})
CREATE (esfbs_bounded_bse:Problem {id:"esfbs_bounded_bse",
  name:"ESFBS with bounded basis-set error",
  definition:"as ESFBS but basis promised to have 1/poly basis-set error for the given V",
  note:"OPEN — posed by arXiv:2103.08215"})
CREATE (esfbs_param_basis:Problem {id:"esfbs_param_basis",
  name:"Electronic structure in parameterized basis",
  definition:"prover supplies both the basis (from a parameterized family) and the witness state",
  note:"OPEN — posed by arXiv:2103.08215"})
CREATE (esfbs_fixed_spin:Problem {id:"esfbs_fixed_spin",
  name:"ESFBS at fixed total spin",
  note:"OPEN — posed by arXiv:2103.08215; the QMA-hardness construction critically uses arbitrary total spin"})
CREATE (mol_gsee:Problem {id:"mol_gsee", name:"Molecular Ground-State Energy Estimation",
  definition:"ground-state energy of the Coulomb Hamiltonian with external potential arising SOLELY from fixed positive point nuclei; freedom = nuclear geometry",
  status:"OPEN — neither QMA-hardness nor BQP membership proven",
  frontier_check:"arXiv sweep 2026-07-16: no resolution found"})
CREATE (nrep_2rdm:Problem {id:"nrep_2rdm", name:"N-Representability of the 2-RDM",
  definition:"does a given 2-RDM extend to a valid N-fermion state (QMA-complete under Turing reductions; Karp-reduction hardness open per 2103.08215 discussion)"})
CREATE (hf_opt:Problem {id:"hf_opt", name:"Lowest-Energy Slater Determinant (Hartree-Fock)",
  definition:"minimize <Phi|H|Phi> over Slater determinants, 1/poly precision"})
CREATE (fhk_eval:Problem {id:"fhk_eval", name:"Hohenberg-Kohn Functional Evaluation",
  definition:"evaluate the universal HK density functional F[rho] (hardness under Turing reductions)"})
CREATE (cure_sign:Problem {id:"cure_sign", name:"Sign-Curing",
  definition:"does a local basis change exist rendering H stoquastic"})
CREATE (chem_sampling:Problem {id:"chem_sampling",
  name:"Classical sampling of quantum-chemistry circuits",
  definition:"sample output distributions of circuit families used in near-term chemical ground-state estimation"})

CREATE (P:ComplexityClass {id:"P", name:"P"})
CREATE (BPP:ComplexityClass {id:"BPP", name:"BPP"})
CREATE (BQP:ComplexityClass {id:"BQP", name:"BQP"})
CREATE (NP:ComplexityClass {id:"NP", name:"NP"})
CREATE (MA:ComplexityClass {id:"MA", name:"MA"})
CREATE (StoqMA:ComplexityClass {id:"StoqMA", name:"StoqMA"})
CREATE (QMA:ComplexityClass {id:"QMA", name:"QMA"})
CREATE (QMA_EXP:ComplexityClass {id:"QMA_EXP", name:"QMA_EXP"})
CREATE (PP:ComplexityClass {id:"PP", name:"PP"})
CREATE (PSPACE:ComplexityClass {id:"PSPACE", name:"PSPACE"})
CREATE (SharpP:ComplexityClass {id:"SharpP", name:"#P"})


// ---------- C. THE VERIFIED REDUCTION CHAIN + HARDNESS LATTICE ----------
CREATE (hk_theorem:Theorem {id:"hk_theorem", name:"Hohenberg-Kohn Theorem (1964)"})

// Root and classic lattice
CREATE (lh5)-[:COMPLETE_FOR {status:"THM", ref:"Kitaev-Shen-Vyalyi 2002 (book)",
  note:"clock construction; QMA verification circuits Karp-reduce to LH5"}]->(QMA)
CREATE (lh2)-[:COMPLETE_FOR {status:"THM", ref:"Kempe-Kitaev-Regev, arXiv:quant-ph/0406180",
  note:"perturbation gadgets"}]->(QMA)
CREATE (lh2_2d)-[:COMPLETE_FOR {status:"THM", ref:"Oliveira-Terhal, arXiv:quant-ph/0504050"}]->(QMA)
CREATE (lh_1d_d12)-[:COMPLETE_FOR {status:"THM",
  ref:"Aharonov-Gottesman-Irani-Kempe, arXiv:0705.4077"}]->(QMA)
CREATE (lh_1d_ti)-[:COMPLETE_FOR {status:"THM", ref:"Gottesman-Irani, arXiv:0905.2419"}]->(QMA_EXP)

// The chain actually used to reach chemistry (verified in 2103.08215)
CREATE (heisenberg_afm)-[:COMPLETE_FOR {status:"THM",
  ref:"Piddock-Montanaro, arXiv:1506.04014 (QIC 2017)",
  note:"holds with poly-bounded nonnegative coefficients", verified:"2026-07-16"}]->(QMA)
CREATE (heisenberg_afm)-[:REDUCES_TO {status:"THM",
  ref:"O'Gorman-Irani-Whitfield-Fefferman, arXiv:2103.08215, Thm 3",
  note:"second-order perturbation theory: large uniform onsite repulsion u0 >= n^(14+3p+2q); t_ij = sqrt(u0*kappa_ij/2)",
  verified:"2026-07-16"}]->(fermi_hubbard)
CREATE (fermi_hubbard)-[:COMPLETE_FOR {status:"THM", ref:"arXiv:2103.08215, Thm 3",
  note:"even with all hopping coefficients of the same sign; no magnetic field (cf. Schuch-Verstraete, who needed a site-specific field)",
  verified:"2026-07-16"}]->(QMA)
CREATE (fermi_hubbard)-[:REDUCES_TO {status:"THM", ref:"arXiv:2103.08215, Sec 4",
  note:"Hubbard instance approximated by ESFBS: orbitals = superpositions of Gaussians; edge {i,j} encoded as a close Gaussian pair at tunable distance gamma_ij; kinetic-energy overlap realizes t_ij; tight beta-Gaussian realizes onsite repulsion; approximation chain H_ES -> H_round -> H_main ~ H_Hubb",
  verified:"2026-07-16"}]->(es_fixedbasis)

// Chemistry-problem statements
CREATE (lh2)-[:REDUCES_TO {status:"THM", ref:"Schuch-Verstraete, arXiv:0712.0483 (Nat Phys 2009)",
  note:"local Hamiltonian encoded in a site-specific external magnetic field; good basis",
  verified:"2026-07-16"}]->(es_v)
CREATE (es_v)-[:COMPLETE_FOR {status:"THM", ref:"arXiv:0712.0483",
  note:"membership: witness = ground state, verifier = QPE; needs the 1/poly promise gap"}]->(QMA)
CREATE (es_fixedbasis)-[:COMPLETE_FOR {status:"THM",
  ref:"O'Gorman-Irani-Whitfield-Fefferman, arXiv:2103.08215, Thm 1 (PRX Quantum 3, 020322)",
  note:"V = 0; instance encoded in the basis choice; incomparable to Schuch-Verstraete (they: magnetic field + good basis; here: no field + artificial basis)",
  verified:"2026-07-16"}]->(QMA)
CREATE (nrep_2rdm)-[:COMPLETE_FOR {status:"THM",
  ref:"Liu-Christandl-Verstraete, arXiv:quant-ph/0609125 (PRL 2007)",
  note:"under Turing reductions"}]->(QMA)
CREATE (hf_opt)-[:COMPLETE_FOR {status:"THM", ref:"arXiv:2103.08215, Thm 2",
  note:"NP-complete for ELECTRONIC-STRUCTURE Hamiltonians specifically (reduction: parameters set so H becomes ~diagonal; encodes independent set); Schuch-Verstraete arXiv:0712.0483 showed NP-hardness for generic quartic fermionic H",
  verified:"2026-07-16"}]->(NP)
CREATE (fhk_eval)-[:HARD_FOR {status:"THM",
  ref:"Schuch-Verstraete arXiv:0712.0483; Whitfield-Schuch-Verstraete 2014",
  note:"under Turing reductions"}]->(QMA)
CREATE (hk_theorem)-[:GUARANTEES_EXISTENCE_OF {status:"THM"}]->(fhk_eval)
CREATE (hk_theorem)-[:DOES_NOT_GUARANTEE_COMPUTABILITY_OF {status:"FACT"}]->(fhk_eval)
CREATE (lh_constgap)-[:COMPLETE_FOR {status:"CONJ",
  note:"quantum-PCP conjecture, open"}]->(QMA)

// Open variants posed by the source paper itself (citable OPEN status)
CREATE (es_fixedbasis)-[:POSES_AS_OPEN {status:"OPEN", ref:"arXiv:2103.08215, Sec 2.1",
  verified:"2026-07-16"}]->(esfbs_bounded_bse)
CREATE (es_fixedbasis)-[:POSES_AS_OPEN {status:"OPEN", ref:"arXiv:2103.08215, Sec 2.1",
  verified:"2026-07-16"}]->(esfbs_param_basis)
CREATE (es_fixedbasis)-[:POSES_AS_OPEN {status:"OPEN", ref:"arXiv:2103.08215, Sec 2.1",
  verified:"2026-07-16"}]->(esfbs_fixed_spin)
CREATE (es_fixedbasis)-[:POSES_AS_OPEN {status:"OPEN", ref:"arXiv:2103.08215, Sec 2.1",
  note:"'Is it still hard when the external potential arises solely from positively charged nuclei at fixed positions?' — this IS Mol-GSEE; open status is now citable",
  verified:"2026-07-16"}]->(mol_gsee)


// ---------- D. CLASS STRUCTURE ----------
CREATE (P)-[:CONTAINED_IN {status:"THM", strictness:"believed strict, unproven"}]->(BPP)
CREATE (BPP)-[:CONTAINED_IN {status:"THM", strictness:"believed strict, unproven"}]->(BQP)
CREATE (BQP)-[:CONTAINED_IN {status:"THM", strictness:"believed strict, unproven"}]->(QMA)
CREATE (NP)-[:CONTAINED_IN {status:"THM", note:"via MA",
  strictness:"believed strict, unproven"}]->(QMA)
CREATE (MA)-[:CONTAINED_IN {status:"THM",
  ref:"Bravyi-DiVincenzo-Oliveira-Terhal, arXiv:quant-ph/0606140",
  strictness:"believed strict, unproven"}]->(StoqMA)
CREATE (StoqMA)-[:CONTAINED_IN {status:"THM", strictness:"believed strict, unproven"}]->(QMA)
CREATE (QMA)-[:CONTAINED_IN {status:"THM", ref:"Marriott-Watrous, arXiv:cs/0506068",
  strictness:"believed strict, unproven"}]->(PP)
CREATE (PP)-[:CONTAINED_IN {status:"THM", strictness:"believed strict, unproven"}]->(PSPACE)


// ---------- E. SIGN PROBLEM (afflicts STOCHASTIC methods only) ----------
CREATE (antisymmetry:Concept {id:"antisymmetry", name:"Fermionic Antisymmetry"})
CREATE (sign_problem:Obstruction {id:"sign_problem", name:"Fermionic Sign Problem"})
CREATE (nodal_bias:Error {id:"nodal_bias", name:"Uncontrolled Nodal Bias"})
CREATE (stoquastic_h:Concept {id:"stoquastic_h", name:"Stoquastic Hamiltonian"})
CREATE (fixednode:Method {id:"fixednode", name:"Fixed-Node Approximation"})
CREATE (fci:Method {id:"fci", name:"Full Configuration Interaction",
  character:"deterministic, exact within basis"})
CREATE (qmc:Method {id:"qmc", name:"Projector / Auxiliary-Field QMC", character:"stochastic"})

CREATE (antisymmetry)-[:INDUCES {status:"FACT", note:"in generic bases"}]->(sign_problem)
CREATE (sign_problem)-[:CAUSES_VARIANCE_BLOWUP_IN {status:"FACT",
  formula:"variance ~ exp(beta*N*Delta_f)"}]->(qmc)
CREATE (sign_problem)-[:DOES_NOT_AFFLICT {status:"FACT"}]->(fci)
CREATE (cure_sign)-[:HARD_FOR {status:"THM",
  ref:"Troyer-Wiese arXiv:cond-mat/0408370; Marvian-Lidar-Hen arXiv:1802.03408; Klassen et al. arXiv:1906.08800"}]->(NP)
CREATE (stoquastic_h)-[:ADMITS {status:"THM", note:"sign-free QMC"}]->(qmc)
CREATE (lh_stoq)-[:COMPLETE_FOR {status:"THM", ref:"arXiv:quant-ph/0606140"}]->(StoqMA)
CREATE (lh_stoq_ff)-[:COMPLETE_FOR {status:"THM", ref:"Bravyi-Terhal, arXiv:0806.1746"}]->(MA)
CREATE (fixednode)-[:REMOVES {status:"FACT"}]->(sign_problem)
CREATE (fixednode)-[:INTRODUCES {status:"FACT"}]->(nodal_bias)


// ---------- F. TENSOR NETWORKS ----------
CREATE (gapped_1d:Concept {id:"gapped_1d", name:"1D Gapped Local Hamiltonian"})
CREATE (area_law:Concept {id:"area_law", name:"Entanglement Area Law"})
CREATE (mps_approx:Concept {id:"mps_approx", name:"Poly-Bond-Dimension MPS Approximation"})
CREATE (mps_global_opt:Problem {id:"mps_global_opt", name:"Global MPS Energy Optimization"})
CREATE (dmrg:Method {id:"dmrg", name:"DMRG"})
CREATE (peps:Method {id:"peps", name:"Projected Entangled Pair States"})

CREATE (gapped_1d)-[:OBEYS {status:"THM", ref:"Hastings, arXiv:0705.2024"}]->(area_law)
CREATE (area_law)-[:IMPLIES {status:"THM"}]->(mps_approx)
CREATE (gapped_1d)-[:GROUND_STATE_COMPUTABLE_IN {status:"THM",
  ref:"Landau-Vazirani-Vidick, arXiv:1307.5143"}]->(P)
CREATE (dmrg)-[:HEURISTICALLY_OPTIMIZES {status:"EMPIRICAL",
  note:"no global optimality guarantee"}]->(mps_approx)
CREATE (mps_global_opt)-[:HARD_FOR {status:"THM", ref:"Eisert, arXiv:quant-ph/0609051"}]->(NP)
CREATE (peps)-[:EXACT_CONTRACTION_HARD_FOR {status:"THM",
  ref:"Schuch-Wolf-Verstraete-Cirac, arXiv:quant-ph/0611050",
  note:"obstructs exact 2D PEPS energy evaluation"}]->(SharpP)


// ---------- G. QUANTUM ALGORITHMS & GUIDED-LH (updated to 2025) ----------
CREATE (qpe:Method {id:"qpe", name:"Quantum Phase Estimation"})
CREATE (overlap_gamma:Concept {id:"overlap_gamma", name:"Trial-State Overlap gamma"})
CREATE (dequantizable:Concept {id:"dequantizable", name:"Classically Dequantizable Regime"})
CREATE (orth_catastrophe:Obstruction {id:"orth_catastrophe", name:"Orthogonality Catastrophe"})
CREATE (adiabatic:Method {id:"adiabatic", name:"Adiabatic State Preparation"})
CREATE (path_gap:Concept {id:"path_gap", name:"Adiabatic Path Spectral Gap",
  note:"can close exponentially fast"})

CREATE (qpe)-[:ESTIMATES {status:"THM", note:"E0 to eps in time poly(1/eps, 1/gamma)"}]->(es_v)
CREATE (qpe)-[:REQUIRES {status:"DEF", note:"efficiency assumes gamma >= 1/poly"}]->(overlap_gamma)
CREATE (guidedlh_poly)-[:COMPLETE_FOR {status:"THM",
  ref:"Gharibian-Le Gall arXiv:2111.09079 (STOC 2022); improved: Cade-Folkertsma-Niesen-Weggemans arXiv:2207.10097; Gharibian-Hayakawa-Le Gall-Morimae arXiv:2207.10250",
  note:"BQP-hardness pushed to 2-local, physically-relevant families and excited states",
  verified:"2026-07-16"}]->(BQP)
CREATE (guidedlh_const)-[:MEMBER_OF {status:"THM", ref:"arXiv:2111.09079",
  note:"constant-precision regime dequantized via classical QSVT sampling techniques",
  verified:"2026-07-16"}]->(dequantizable)
CREATE (guidedlh_stoq)-[:HARD_FOR {status:"THM", ref:"arXiv:2509.25829 (2025)",
  note:"guided LH restricted to stoquastic H is (promise) BPP-hard — classically as hard as it gets short of full BQP",
  verified:"2026-07-16"}]->(BPP)
CREATE (guidable_lh)-[:STUDIED_IN {status:"THM",
  ref:"Weggemans-Folkertsma-Cade arXiv:2302.11578",
  note:"Merlinized guidable variants; implications for heuristic ansatz preparation and quantum PCP",
  verified:"2026-07-16"}]->(lh_constgap)
CREATE (dmrg)-[:PROVIDES_GUIDING_STATES_FOR {status:"EMPIRICAL",
  ref:"physically-motivated guiding states: arXiv:2509.25815; classical/dequantized algorithms: arXiv:2411.16163, arXiv:2409.04161",
  verified:"2026-07-16"}]->(guidedlh_poly)
CREATE (orth_catastrophe)-[:DRIVES {status:"EMPIRICAL",
  note:"gamma shrinks toward exp(-cN) generically at large N; central contested premise of the advantage debate"}]->(overlap_gamma)
CREATE (overlap_gamma)-[:VANISHING_VOIDS {status:"COROLLARY"}]->(guidedlh_poly)
CREATE (adiabatic)-[:ASSUMES {status:"DEF", note:"1/poly minimum gap along path"}]->(path_gap)


// ---------- G2. SAMPLING-HARDNESS AXIS (new, 2025) ----------
CREATE (chem_sampling)-[:CONDITIONALLY_HARD {status:"THM",
  ref:"arXiv:2503.21041 (Berkeley 2025); arXiv:2504.12893 (Tokyo/QunaSys 2025)",
  note:"classical simulation hardness for chemistry-circuit families under generalized complexity conjectures; NOTE: a SAMPLING statement about circuits, distinct from decision-problem (energy-estimation) hardness — does not resolve mol_gsee",
  verified:"2026-07-16"}]->(SharpP)


// ---------- H. TRACTABLE ISLANDS & THE EVIDENCE DEBATE (updated to 2026) ----------
CREATE (quadratic_h:Concept {id:"quadratic_h", name:"Quadratic Fermionic Hamiltonian"})
CREATE (worst_case_caveat:MetaNote {id:"worst_case_caveat",
  name:"QMA-hardness is a worst-case statement over engineered instances",
  support:"arXiv:2103.08215 states its basis is 'artificial' and that theory/numerics suggest physically realistic potentials always admit a good poly-size basis"})
CREATE (physical_molecules:Concept {id:"physical_molecules", name:"Physically Occurring Molecules"})
CREATE (advantage_evidence:Concept {id:"advantage_evidence",
  name:"Evidence debate: exponential quantum advantage in ground-state chemistry",
  key_refs:"Lee et al. arXiv:2208.02199 (Nat Commun 2023, skeptical); FeMoco model solved classically to chemical accuracy arXiv:2601.04621 (Jan 2026); Fe4S4 CAS(54,36) Blackwell-DMRG benchmark arXiv:2603.28648 (Mar 2026)"})
CREATE (avg_case_tractability:Concept {id:"avg_case_tractability",
  name:"Average-Case Tractability of Molecular Instances", note:"unproven in either direction"})

CREATE (quadratic_h)-[:MEMBER_OF {status:"THM",
  ref:"Valiant matchgates arXiv-era 2002; Terhal-DiVincenzo arXiv:quant-ph/0108010"}]->(P)
CREATE (es_v)-[:HARDNESS_IS {status:"FACT"}]->(worst_case_caveat)
CREATE (es_fixedbasis)-[:HARDNESS_IS {status:"FACT", verified:"2026-07-16"}]->(worst_case_caveat)
CREATE (physical_molecules)-[:NOT_PROVEN_HARD_INSTANCE_OF {status:"FACT"}]->(es_fixedbasis)
CREATE (advantage_evidence)-[:STRENGTHENS_CLASSICAL_SIDE_VIA {status:"EMPIRICAL",
  note:"canonical advantage targets (FeMoco, Fe4S4) now solved/benchmarked classically to chemical accuracy; DMRG posited as mandatory classical reference for any advantage claim",
  verified:"2026-07-16"}]->(avg_case_tractability)
CREATE (dmrg)-[:SUPPLIES_BENCHMARKS_FOR {status:"EMPIRICAL",
  ref:"arXiv:2603.28648", verified:"2026-07-16"}]->(advantage_evidence)


// ---------- I. METHODS vs CERTIFICATION ----------
CREATE (krylov:Method {id:"krylov", name:"Krylov / Lanczos Iteration"})
CREATE (rayleigh_ritz:Method {id:"rayleigh_ritz", name:"Rayleigh-Ritz on Krylov Subspace"})
CREATE (temple_lehmann:Method {id:"temple_lehmann", name:"Temple / Lehmann-type Bounds"})
CREATE (two_sided_cert:DesignPrinciple {id:"two_sided_cert",
  name:"Two-Sided Certification (krylov-quantum-solver)"})

CREATE (fci)-[:SOLVES_EXACTLY_IN_BASIS {status:"DEF"}]->(es_fixedbasis)
CREATE (fci)-[:COST_SCALES_AS {status:"FACT", note:"poly(dim Fock) = exp(N)"}]->(fock)
CREATE (krylov)-[:MATVEC_COST_SCALES_AS {status:"FACT",
  note:"O(dim Fock * poly(N)) per matvec"}]->(fock)
CREATE (rayleigh_ritz)-[:CERTIFIES_UPPER_BOUND_ON {status:"THM", ref:"variational principle",
  note:"unconditional"}]->(es_fixedbasis)
CREATE (temple_lehmann)-[:CERTIFIES_LOWER_BOUND_ON {status:"THM",
  note:"conditional on a gap estimate as input"}]->(es_fixedbasis)
CREATE (two_sided_cert)-[:IS_RATIONAL_RESPONSE_TO {status:"DESIGN",
  note:"worst-case QMA-completeness in the solver's exact regime (fixed basis) + unproven instance-hardness + intensifying classical-heuristic frontier -> verify, don't trust"}]->(worst_case_caveat)


// ---------- J. META / FLAGS ----------
CREATE (amorphous_bc:MetaNote {id:"amorphous_bc", name:"'Amorphous boundary conditions'",
  note:"no established complexity-theoretic result attaches to this term"})
CREATE (translational_invariance:Concept {id:"translational_invariance", name:"Translational Invariance"})
CREATE (this_graph:MetaNote {id:"this_graph", name:"This Knowledge Graph",
  note:"models argument structure; is not itself a proof; v2 arXiv-verified 2026-07-16"})
CREATE (false_claim:MetaNote {id:"false_claim",
  name:"'The exact complexity of Mol-GSEE is known'",
  note:"FALSE — open status citable via arXiv:2103.08215; frontier sweep 2026-07-16 found no resolution"})

CREATE (translational_invariance)-[:NEAREST_RIGOROUS_AXIS_TO {status:"META",
  note:"TI-1D LH is QMA_EXP-complete (lh_1d_ti)"}]->(amorphous_bc)
CREATE (this_graph)-[:FLAGS {status:"META"}]->(false_claim)


// ============================================================
// EXAMPLE QUERIES — run separately AFTER loading.
// ============================================================
// The verified reduction chain, root to chemistry:
//   MATCH p = (a)-[:REDUCES_TO*1..4]->(b {id:"es_fixedbasis"}) RETURN p
// Everything that is not a proven theorem:
//   MATCH ()-[r]->() WHERE r.status <> "THM" RETURN r
// All session-verified edges:
//   MATCH ()-[r]->() WHERE r.verified IS NOT NULL RETURN r
// All open problems and who posed them:
//   MATCH (a)-[r:POSES_AS_OPEN]->(b) RETURN a.name, b.name, r.ref
// The evidence-debate subgraph:
//   MATCH (n {id:"advantage_evidence"})-[r]-(x) RETURN n, r, x
