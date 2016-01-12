# kidney\_solver

kidney\_solver is a program for the _kidney-exchange problem_, using Python 2 and the Gurobi IP solver.

## Prerequisites

- Python 2
- [Gurobi](http://www.gurobi.com)

## Usage

The program reads an instance from standard input, and writes the result to standard output. If no non-directed donors (NDDs) are used, the input format is the `.input` format generated by, for example, David Abraham's Saidman generator. The example in the `example_data` directory has 64 donor-patient pairs numbered 0, ..., 63, and has 1025 edges. The rows representing edges are in _source target weight_ format. The final row of the input contains -1 three times.

If NDDs are used, the input described above is followed by additional data in a very similar format. The `.ndds` file in the `example_data` directory, for example, has 6 NDDs and 188 edges from NDDs. NDDs are numbered from zero. The edges are in _source-NDD target-pair weight_ format.

The program `utils/convert.py` can be used to convert from `.wmd` format.

The `kidney_solver` program has three required command-line arguments: cycle cap, chain cap, and formulation. Note that the chain cap is the maximum permitted number of edges in a chain, _excluding the dummy arc to the NDD_. The formulation can be:

- uef: Edge formulation with unrestricted cycle and chain sizes
- hpief\_prime: A slightly-modified hybrid PIEF
- picef: Position-indexed chain-edge formulation
- cf: Cycle formulation, with one variable per cycle or chain

In addition, the optional flag `-r` can be used to solve on a copy of the graph with vertices relabelled in descending order of out-degree plus in-degree, which may result in a smaller IP model. To set a time limit of LIMIT seconds, use `-t LIMIT`.

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
