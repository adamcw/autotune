# MANUAL #

## INTRODUCTION ##

The library `qc` is designed to simulate and be used to benchmark a quantum
computer running error correction that makes use of minimum weight perfect
matching. This includes the surface code and 3-D topological cluster states. It
allows the user to simulate any error model that can be described by a single
integer per qubit, with the caveat that only single- and two-qubit physical
gates are implemented in the current version. A flexible Pauli error based
library `depolar` is provided which should be sufficient for most purposes. The
current version is well-suited to obtaining architecture specific threshold
error rates and the performance of logical memory constructed using both
boundaries and defects. Later versions will include the ability to benchmark
braided logical gates.

It is intended that the user will use qc by rewriting `ex1.c` in `depolar` in a
manner tailored to their specific architecture, set of error models and
problem. This document describes, step-by-step, the code that needs to be
written. 

This distribution is coupled to the `Blossom V` library, implementing the
minimum-weight perfect matching algorithm. Wrapping this library is `match`
which acts as a translation layer between the `Autotune` library and `Blossom
V`.

## WHAT'S IN AUTOTUNE? ##

`Autotune` is comprised of three modules. 

- The `Autotune` libraries known as `libautotune`. This contains the following libraries:
    - `match` (A wrapper for the Blossom V matching library)
    - `depolar` (Depolarizing Quantum Computer)
    - `qc` (Quantum Computer)
    - `bfs` (Breadth First Search)
    - `bheap` (Binary Heap)
    - `cdllist` (Circular Double-Linked List)
    - `dllist` (Double-Linked List)
    - `hash_table` (A Hash Table)
    - `llist` (Linked List)
    - `memory` (Memory functions)
    - `my_time` (Timing functions)
    - `random` (A Mersenne Prime Twister random number generator library)
- Examples of the use of the `Autotune` library This contains the following examples:
    - `ex1` (The use of `Autotune` to implement the Surface Code on a
        solid-state quantum computer)
- Processing scripts for generating graphs based on the output of `Autotune`

`Autotune` additionally makes use of an external library `Blossom V` which must
be downloaded seperately and installed. See `README.md` for more details.

## GATES AND ERROR MODELS ##

`depolar` provides implementations of a range of standard gates:

    void dp_init_X(DP_QC *dp_qc, QUBIT *q);
    void dp_init_Z(DP_QC *dp_qc, QUBIT *q);
    void dp_iden(DP_QC *dp_qc, QUBIT *q, int gate);
    void dp_iden_init_X(DP_QC *dp_qc, QUBIT *q);
    void dp_iden_init_Z(DP_QC *dp_qc, QUBIT *q);
    void dp_iden_H(DP_QC *dp_qc, QUBIT *q);
    void dp_iden_cnot(DP_QC *dp_qc, QUBIT *q);
    void dp_iden_cZ(DP_QC *dp_qc, QUBIT *q);
    void dp_iden_swap(DP_QC *dp_qc, QUBIT *q);
    void dp_iden_meas_X(DP_QC *dp_qc, QUBIT *q);
    void dp_iden_meas_Z(DP_QC *dp_qc, QUBIT *q);
    void dp_iden_transform_e(QUBIT *q);
    void dp_dead(DP_QC *dp_qc, QUBIT *q, int gate);
    void dp_dead_init_X(DP_QC *dp_qc, QUBIT *q);
    void dp_dead_init_Z(DP_QC *dp_qc, QUBIT *q);
    void dp_dead_H(DP_QC *dp_qc, QUBIT *q);
    void dp_dead_cnot(DP_QC *dp_qc, QUBIT *q);
    void dp_dead_cZ(DP_QC *dp_qc, QUBIT *q);
    void dp_dead_swap(DP_QC *dp_qc, QUBIT *q);
    void dp_dead_meas_X(DP_QC *dp_qc, QUBIT *q);
    void dp_dead_meas_Z(DP_QC *dp_qc, QUBIT *q);
    void dp_H(DP_QC *dp_qc, QUBIT *q);
    void dp_cnot(DP_QC *dp_qc, QUBIT *q1, QUBIT *q2);
    void dp_cZ(DP_QC *dp_qc, QUBIT *q1, QUBIT *q2);
    void dp_swap(DP_QC *dp_qc, QUBIT *q1, QUBIT *q2);
    int dp_meas_X(DP_QC *dp_qc, QUBIT *q, SET *set1, SET *set2);
    int dp_meas_Z(DP_QC *dp_qc, QUBIT *q, SET *set1, SET *set2);

For the most part, these functions do what would naively be expected from their
names -- given a `dp_qc` (depolarizing quantum computer) and a qubit or qubits
within that computer, they apply the named gate, including stochastic errors.
The dead gate should only be used after a measurement and before an
initialization of a qubit. The stochastic error models are described in files
`gate_name_em`. The format of these files is:

* number of qubits the gate is applied to (1 or 2) scale factor controlling the
* relative strength of different gate
error models (arbitrary float)
* number of lines in the error model a list of relative error strengths and
* error values (two or three
integers)
* the gate duration (integer)

For example, `cnot_em` contains:

    2
    1.0
    15
    1 0 1
    1 0 2
    1 0 3
    1 1 0
    1 1 1
    1 1 2
    1 1 3
    1 2 0
    1 2 1
    1 2 2
    1 2 3
    1 3 0
    1 3 1
    1 3 2
    1 3 3
    1

In `depolar`, the error models all involve only Pauli operators, with 0=I, 1=X,
2=Z, 3=Y. The fourth line above therefore represents an IX error with relative
strength 1. Since the relative strengths sum to 15 and are all equal, this is
the standard depolarizing channel. If one wishes to strengthen or weaken all
cnot errors, this can be done by increasing or decreasing the scale factor.
`depolar` should be consulted if the user wishes to create their own gates or
error models. H is shipped with no errors and no duration, along with a
function `iden_H` that also has no errors and no duration, to allow one to
implement `meas_X` in a more physical H then `meas_Z` manner without
introducing extra errors. This can easily be modified by using the X, Y or Z
error models.


## SETS ##

The measurement functions listed above take two set arguments. In the surface
code example, a set is usually a pair of consecutive stabilizer measurements.
Each stabilizer measurement in the surface code has value +1 or -1. If the
product of measurement values associated with a set is -1, an error must have
occurred locally.  Before performing measurements, the user must create sets
using:

    SET *qc_create_set(int type, int i, int j, int num_meas_left, SET *bdy);

`type` is either `PRIMAL` or `DUAL`. In the surface code example, `PRIMAL` set
are associated with X-stabilizer measurements and `DUAL` with Z-stabilizer
measurements.

`i`, `j` are for debugging only -- it makes sense to set these to values
representing the physical location of the stabilizer measurement, but the user
could equally well set them to arbitrary values.

`num_meas_left` is the number of measurements left which in the surface code
example is always 2 since every set contains two measurements. Note that the
set does not know which measurements will be associated with it, only how
many.

`bdy` allows one to specify a nearby boundary set. Sets on the boundary of the
lattice need to be told which boundary they are near to ensure that if the set
detects an error, it is possible to match that error to the appropriate
boundary. Sets not on the boundary of the lattice should have this argument set
to NULL.

`ex1.c` in `depolar` contains code showing how to set up appropriate boundary
sets (in `create_sc_dp_qc()`: create surface code depolarizing quantum
computer) and continue creating appropriate sets as the computation proceeds
(in `measure_stabilizers()`).


## SYNDROMES ##

A syndrome is created for every location a primal or dual stabilizer
measurement can be performed:

    SYNDROME *qc_create_syndrome();

It may seem surprising that no arguments are required. The purpose of the
syndrome structure is to keep track of the progress of the error correction.
Syndromes are associated with sets when they are close to completion
(`num_meas_left` close to 0):

    void qc_associate_syndrome(SET *set, SYNDROME *syn);

and unassociated when the set is complete:

    void qc_unassociate_syndrome(SYNDROME *syn);

Whenever every syndrome has been associated with a completed set, a round of
the error correction has been completed and a variable `big_t` is incremented.
`create_sc_dp_qc()` and `measure_stabilizers()` show how to use syndromes.


## NESTS AND LATTICES ##

With code as described above written, `qc` is capable of constructing a graph
which we call a nest of balls and sticks with each ball corresponding to a set
and each stick representing the total probability that a single error leads to
a pair of detection events in the linked sets/balls. All of this is internal.
Before matching can occur, the nest needs to be converted into a leaner graph
which we call a lattice of dots and lines using:

    void qc_convert_nests(QC *qc, int undo);

The undo flag (TRUE/FALSE) controls whether information is recorded to allow
the conversion to be undone. Vertices are created and associated with dots that
are ultimately associated with stochastic detection events. The function:

    void qc_trim_nests(QC *qc, int big_t);

Is used to discard old nest data to ensure that the nest memory requirements
remain finite.


## MATCHING ##

Given a lattice, the weight of a line is `-ln(p_line)`. The weight is a positive
number that grows larger as `p_line` grows smaller. Given vertices associated
with a selection of dots, the weight of an edge connecting a pair of vertices
is the minimum weight path connecting the associated pair of dots. Given a
lattice with associated vertices, our full minimum weight perfect matching
function returns a list of pairs of vertices such that every vertex is paired
with exactly one other vertex, and the total weight of all pairing edges is
minimal. The `qc` function:

    void qc_mwpm(QC *qc, int undo);

Calls the `match` function:

    void m_mwpm(MATCHING *m, int undo);

which simply chooses a small number of edges at random. This is sufficient to
test that code is working. The edges are returned as a list of structures
internal to the matching of type `AUG_EDGE`. The `ex1.c` function:

    void correct_mts(SC_DP_QC *sc_dp_qc);

Provides an example of how to use these augmented edges to apply corrections to
the measurement data. Essentially, a Pauli frame is updated and a path of
measurements between the vertices is flipped.  Furthermore, lines of
measurements starting at the vertices and leading into the future are flipped,
creating a U-shape of flipped measurements. This is the correct thing to do for
the surface code.  For 3-D topological cluster states, one would only flip the
path of measurements between the vertices.

Logical errors are detected using:

    void test_correct(SC_DP_QC *sc_dp_qc);

The current `sc_dp_qc` is copied, two more rounds of stabilizers are measured
on this copy, the copied nest is converted into a lattice recording information
to allow this operation to be undone (matchings are not copied), matching is
performed recording information to undo, corrections are applied to the copy of
`sc_dp_qc`, and the errors on the qubits and Pauli frames are then examined to
see if a logical error has occurred since the last `test_correct()` call.


## ABOUT T_DELETE ##

`t_delete` is an essential parameter to the running of Autotune. It tells
Autotune how many timesteps to keep in memory for the sake of error correction.
High error rates (-p) require more memory to be kept to retain accuracy. As a
rough guide we recommend a `t_delete` of around 100 to be safe for values of p up
to around 0.05. `t_delete` can rapidly need to be very large as the error rate
approaches the threshold. At lower error rates, `t_delete` can be 10 or even
lower with little impact on the measured logical error rate. 

The use of a `t_delete` that is too small for a given problem can lead to a
slightly elevated logical error rate. It is therefore recommended, where
accuracy is paramount, to test at a given error rate for a number of different
values of `t_delete` and determine which is best suited.

For more information on choosing a `t_delete`, please contact us.

## RUN PROCESSING ##

See process/README.md
