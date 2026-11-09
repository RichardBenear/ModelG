# MG12-H0 initial specification and static checks — MG105

We have made the proposed hypergraph translation precise enough to build: particles and resources begin as named token states attached to relational sites, reactions replace local state relations, and the new tick model prevents every participant—including unchanged catalysts—from acting again within that tick. The local reaction schemas preserve the existing finite-count rates and token accounting, and their behavior is independent of arbitrary token labels. This is a specification and static-check result, not a working hypergraph dynamics engine or evidence of particle persistence, emergent dimension or a quantum substrate. Transport execution, stochastic evolution and numerical agreement remain to be implemented and tested.

Richard approved all five MG104 direction choices: MG12 proceeds while G5 remains closed/unresolved and its larger follow-up is deferred; named species remain explicit initially; catalysts obey one-use-per-tick; dimension is investigated relationally; finite resources precede regeneration/quantum extensions.

The versioned contract is spec/MG12H0_contract.md, published at3f37eab before the static checks. It specifies typed Site/Token incidence, species and attachment replacement, fixed-reference A/B rates, spent Z/W retention, pairwise transport semantics and two separate executors: the existing-law continuous-time reference and a proposed strict snapshot-tick matching process.

## Completed static evidence

- 101exact rational/count-increment assertions passed on four tiny reaction fixtures, with at most seven tokens.
- 5,164exhaustive bijective relabelings preserved match roles, rates and output effects after mapping labels back.
- The two X catalysts in R4 are distinct unordered participants; per-match rate2/Omega_count² produces the required total nX(nX−1)nY/Omega_count², rather than half or twice that rate.
- All local reaction increments preserve seven-species total count.
- Participant-disjoint reaction replacements commute on the fixtures.
- An unchanged R4 catalyst overlaps and therefore conflicts with an X-loss event involving that token.
- X created by G→X is absent from the original snapshot's feedback matches, even though it can participate in a later snapshot.

Read contract_checks.json and src/check_mg12h0_contract.py for exact fixtures/results. The separate arithmetic path is an independently written table of count propensities and integer species increments within the same script. It shares the contract and author; it is not independent engine replication. Source/contract SHA256 values refer to local UTF8 bytes; Git line-ending normalization may require reproducing those original bytes.

## What is deliberately not claimed

The static matcher represents a single site's named token states; the full incidence-hypergraph storage and matcher do not yet exist. The proposed transport executor, link-sharing behavior, event clocks, rejected proposals, multi-tick ancestry and finite-delta statistics have not been tested. No random draws, time integrations or pattern experiments were run. These checks do not constitute full H1 qualification or an emergent substrate result.

The strict tick executor samples frozen per-match exponential proposals, selects participant-disjoint proposals in time order and commits together. It is a specified conflict model, not a uniform draw over compatible event sets. It has a first-order correspondence to the old generator for fixed finite states; finite-tick kinetics can differ, especially because catalyst reuse is now restricted. Those differences must be measured rather than hidden by referring to the old kernel as equivalent.

## Next implementation step

Build the tiny explicit-incidence reference engine and its exact/strict-tick executors, then qualify reaction/transport matching, replay, labeling invariance, catalyst saturation and clock refinement against independent count references. Before any numerical development, publish a fixed fixture/seed plan, source hashes, development cap and accuracy criteria. The contract reserves a proposed300seconds single-worker development scope to concretize; no scientific production budget, sample count or full-lattice run is yet registered. Final production approval follows a measured ready-to-run packet. Earlier closed ModelG/CRR results remain unchanged.
