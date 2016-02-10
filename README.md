# kidney\_solver

kidney\_solver is a program for the _kidney-exchange problem_, using Python 2 and the Gurobi IP solver.

## Prerequisites

- Python 2
- [Gurobi](http://www.gurobi.com)
- Nose, for the tests. Run ``nosetests`` from the base directory.


## Sources of IP formulations

- Extended edge formulation: Miguel Constantino, Xenia Klimentova, Ana Viana, Abdur Rais, [New insights on integer-programming models for the kidney exchange problem](http://www.sciencedirect.com/science/article/pii/S0377221713004244), European Journal of Operational Research, Volume 231, Issue 1, 16 November 2013, Pages 57-68.
- PICEF and HPIEF: Forthcoming work by John Dickerson, David Manlove, Benjamin Plaut, Tuomas Sandholm, and James Trimble.
- Hybrid extended edge formulation, using PICEF for chains: Xenia Klimentova.
- Cycle formulation: Roth, Alvin E., Tayfun Sönmez, and M. Utku Ünver. 2007. [Efficient Kidney Exchange: Coincidence of Wants in Markets with Compatibility-Based Preferences.](https://www.aeaweb.org/articles.php?doi=10.1257/aer.97.3.828), American Economic Review, 97(3): 828-851 and David J. Abraham, Avrim Blum, and Tuomas Sandholm, [Clearing Algorithms for Barter Exchange Markets:
Enabling Nationwide Kidney Exchanges](http://www.cs.cmu.edu/~dabraham/papers/abs07.pdf), in Proceedings of ACM-EC 2007: the Eighth ACM Conference on Electronic Commerce.

## Usage

The program reads an instance from standard input, and writes the result to standard output. If no non-directed donors (NDDs) are used, the input format is the `.input` format. The example in the `example_data` directory has 64 donor-patient pairs numbered 0, ..., 63, and has 1025 edges. The rows representing edges are in _source target weight_ format. The final row of the input contains -1 three times.

If NDDs are used, the input described above is followed by additional data in a very similar format. The `.ndds` file in the `example_data` directory, for example, has 6 NDDs and 188 edges from NDDs. NDDs are numbered from zero. The edges are in _source-NDD target-pair weight_ format.

The program `utils/convert.py` can be used to convert from `.wmd` format. (Generated instances in this format can be found on [PrefLib](http://www.preflib.org/data/matching/kidney/).)

The `kidney_solver` program has three required command-line arguments: cycle cap, chain cap, and formulation. Note that the chain cap is the maximum permitted number of edges in a chain, _excluding the dummy arc to the NDD_. The formulation can be:

- ``uef``: Edge formulation with unrestricted cycle and chain sizes
- ``hpief_prime``: A hybrid PIEF with an additional refinement that avoids the need for variables in first position
- ``hpief_prime_full_red``: ``hpief_prime``, with further reduction by generating cycles
- ``hpief_2prime``: ``hpief_prime``, with an additional refinement that avoids the need for variables in position equal to the cycle cap. Note that if the cycle cap is less than 3, ``hpief_prime`` is used instead.
- ``hpief_2prime_full_red``: ``hpief_2prime``, with further reduction by generating cycles
- ``eef``: Reduced extended edge formulation (with a slight modification to the symmetry-breaking constraints)
- ``eef_full_red``: ``eef`` with further reduction by generating cycles
- ``picef``: Position-indexed chain-edge formulation
- ``cf``: Cycle formulation, with one variable per cycle or chain

The optional flag `-r` can be used to solve on a copy of the graph with vertices relabelled in descending order of out-degree plus in-degree, which may result in a smaller IP model. To set a time limit of LIMIT seconds, use `-t LIMIT`.

If the cycle formulation or PICEF is used, failure-aware matching with uniform edge failure probability can be performed with `-p EDGE-SUCCESS-PROB`.

*Example 1:* .wmd format input

```
python utils/convert.py < example_data/MD-00001-00000100.wmd | python kidney_solver/kidney_solver.py 3 3 picef
```

*Example 2:* input in .input and .ndds format

```
cat example_data/MD-00001-00000100.input example_data/MD-00001-00000100.ndds | python kidney_solver/kidney_solver.py 3 3 picef
```

## Output

The output should be mostly self-explanatory. Each row for the `cycles` listing is a list of donor-patient pair indices. Each row of the `chains` listing is the NDD index, followed by a list of donor-patient pair indices.

## Utility to count cycles and chains

A Python utility for counting cycles and chains in an instance is also included. This reads from standard input and takes the cycle and chain caps as command-line arguments. Example usage:

```
cat example_data/MD-00001-00000100.input example_data/MD-00001-00000100.ndds | python kidney_solver/count_cycles_and_chains.py 3 3
```

Note that this will probably run quite a bit quicker if you use Pypy rather than CPython.

## Utility to sparsify instances

The `sparsify.py` instance can be used to delete each edge from an instance with some given probability. The program reads from standard input in the `.input` + `.ndds` format, and writes to standard output in the same format. The probability that each edge will be _kept_ is a command-line argument.

```
cat example_data/MD-00001-00000100.input example_data/MD-00001-00000100.ndds | python kidney_solver/sparsify.py .05
```

## Utilities to randomise edge weights

The `utils/add_random_real_weights.awk` tool sets the weights on edges to random numbers in the range [0,1). Set the random-number seed on the command line using the `seed` variable:

```
cat example_data/MD-00001-00000100.input | awk -f utils/add_random_real_weights.awk -v seed=$RANDOM
```

The `utils/add_random_integer_weights.awk` tool sets the weights on edges to random integers in the range [`lower`, `upper`), where the variables `lower` and `upper` are set on the command line. Set the random-number seed on the command line using the `seed` variable:

```
cat example_data/MD-00001-00000100.input | awk -f utils/add_random_integer_weights.awk -v seed=$RANDOM -v lower=5 -v upper=10
```

## Alternatives

This is a (probably very incomplete) list of other software for kidney exchange.

- https://github.com/rma350/kidneyExchange
- https://github.com/JohnDickerson/KidneyExchange
- https://github.com/ptoulis/kidney-exchange

## Contact

I'd be more than happy to try to answer any questions: james.trimble at yahoo.co.uk

