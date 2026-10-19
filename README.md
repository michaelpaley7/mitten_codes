# Evidence pack: undetectable logical fault sets in the released mitten schedules

Each file in `witnesses/` lists an explicit set of single-qubit Pauli faults for one released schedule and one half
of the syndrome-extraction cycle. Each set fires no detector and flips a logical observable. `check_witness.py`
verifies any of them using **only the authors' own files** (from github.com/a7b/yarn, commit 1b58edd) and Stim.

## Run it (about a minute)

```
git clone https://github.com/a7b/yarn && git -C yarn checkout 1b58edd
pip install numpy stim
python -I check_witness.py "yarn/processor_codes/mitten/[[200,40,12]]" witnesses/mitten_200_40_12_Z.json
```

Expected output, abbreviated:

```
{"faults": 9, ..., "detectors_fired": 0, "logicals_flipped": 23, "result": "PASS: undetectable logical fault set"}
```

Options:
- `--cycles 2 --fault-cycle 0` gives a two-cycle memory.
- `--drop-one` is a negative control: removing any single fault makes detectors fire.

## What the checker does

1. It reads `Hx.npy`, `Hz.npy`, the authors' logical operators `Lx.npy` and `Lz.npy`, and
   `hook_free_SE_cycle_schedule.json`.
2. It confirms that each check's gate sequence in the JSON reproduces its row of Hx or Hz exactly.
3. It confirms that every listed ancilla fault spreads to exactly the data qubits that its check touches later, in
   the JSON order.
4. It builds a memory circuit that executes the JSON layers in the listed order, for several cycles:
   - first the 12 X-check layers (CNOT ancilla→data; ancilla in |+⟩, measured in X);
   - then the 12 Z-check layers (CNOT data→ancilla; ancilla in |0⟩, measured in Z).
5. It confirms that the noiseless circuit is deterministic.
6. It inserts the faults with probability 1 and reports how many detectors fire and how many logicals flip.

For an "X" witness (X faults), the checker uses a Z-basis memory and the authors' `Lz`. For a "Z" witness, it uses an
X-basis memory and `Lx`.

## Results (all 26 pass in 2- and 3-cycle memories; every drop-one control fails)

| Schedule | X half | Z half | Code distance (paper) |
|---|---|---|---|
| mitten [[150,30,10]] | 8 | 8 | 10 (Table VIII reports 8) |
| mitten [[200,40,12]] | **9** | **9** | 12 |
| mitten [[300,60,14]] | **13** | **13** | 14 |
| mitten [[500,100,16]] | 16 | **9** | 16 |
| mitten [[540,108,18]] | 18 | **13** | 18 |
| mitten [[630,126,≤20]] | 14 | 13 | ≤ 20 |
| mitten [[780,156,≤22]] | 17 | 16 | ≤ 22 |
| mitten [[975,195,≤24]] | 23 | 20 | ≤ 24 |
| structured [[300,60,9]] (repository folder) | 9 | 9 | — |
| structured [[330,66,≤12]] | 7 | 7 | ≤ 12 |
| structured [[600,120,≤14]] | 11 | 10 | ≤ 14 |
| structured [[600,120,≤16]] | 10 | 11 | ≤ 16 |
| structured [[840,168,≤18]] | 13 | 13 | ≤ 18 |

These are fault counts of explicit sets, so they are upper bounds on the circuit distance under this circuit
reading, not exact distances. Bold marks values below a code distance that the paper states as exact.
