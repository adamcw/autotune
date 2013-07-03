# README #

## Installing Autotune ##

`Autotune` and its dependent libraries can be downloaded and installed by
calling `make install` from within the `autotune` directory. 

This procedure will download and install `Blossom V` a software library
written by Vladimir Kolmogorov and required by `Autotune`. Please take careful
note of the license of this software as outlined in `blossomv/LICENSE.TXT`
after installation. Research using `Autotune` and therefore `Blossom V` must
cite `Blossom V` in any publications, and commercial use is forbidden without
consent. Please email `vnk@ist.ac.at` with issues regarding licensing of
`Blossom V`.

### Building Autotune ###

To build (make) `Autotune`, `Blossom V`, and the example files:

    cd ./autotune;
    make install;

Now you can run any of the examples from the directory:

    autotune/ex/bin

### Complete Installation Instructions ###

    git clone git@github.com:adamcw/autotune.git
    cd ./autotune/
    make install
    cd ./process
    ./run.py --input example_runs.json

## Issues ##

If there are build issues, you can try and build using a specific version of
`gcc` by modifying your environment to export `$CXX`. If using bash, edit your
`.bashrc` or `.bash_profile` to include:

    export CXX=g++;

## What's in Autotune? ##

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
    - `polyestimate` (A library for estimating logical error rates given a
        database of pre-simulated results)
    - `random` (A Mersenne Prime Twister random number generator library)
- The `tools` that make use of `Autotune`. Exist as advanced examples that
  provide a useful function beyond demonstrating what `Autotune` can do.
    - `polyestimate` (A tool for estimating logical error rates given a
        database of pre-simulated results)
- Examples of the use of the `Autotune` library This contains the following examples:
    - `ex1` (The use of `Autotune` to implement the Surface Code on a
        solid-state quantum computer)
- Processing scripts for generating graphs based on the output of `Autotune`

## Building Autotune ##

You can compile `libautotune` on its own:

    make libautotune;

You can compile `blossomv` on its own:

    make blossomv;

You can compile the examples on their own (though require `libautotune` and
  `blossomv` to be compiled first)

    make ex;

You can compile just the processing pipeline:

    make process;

You can compile `libautotune` and the examples `ex` at the same time (helpful if
  modifying `libautotune` and examples simultaneously through development):

    make autotune;

You can compile all targets at once:

    make all;

## Processing Autotune ##

    See: process/README.md

## Polyestimate Tool ##

    See: tools/polyestimate/README.md

    Austin G. Fowler. "Polyestimate: instantaneous open source surface code analysis"
    Awaiting Publication (2013)
    Note: arXiv:1307.0689, http://topqec.com.au/autotune.html

## More Information ##

    See: MANUAL

## Contact ##

    Austin Fowler (afowler@topqec.com.au)
    Adam Whiteside (awhiteside@topqec.com.au)

## Website ##

    Topological Quantum Error Correction
    http://www.topqec.com.au

## License ##

Copyright 2013 Austin G. Fowler, Adam C. Whiteside, Angus L. McInnes, Alimohammad Rabbani

This software cannot be used for commercial purposes, or distributed without
permission from the author. This software is made available for research use
only. Please contact us for alternate licensing arrangements. 

You can cite `Autotune` using the following citation:

    A. Fowler, A. Whiteside, A. McInnes, and A. Rabbani. 
    "Topological code Autotune." Phys. Rev. X 2, 041003 (2012) 
    Note: arXiv:1202.6111, http://topqec.com.au/autotune.html

You can cite the `Polyestimate` tool using the following citation:

    Austin G. Fowler. "Polyestimate: instantaneous open source surface code analysis"
    Awaiting Publication (2013)
    Note: arXiv:1307.0689, http://topqec.com.au/autotune.html


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
