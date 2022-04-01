"""

MPI python script to calculate DTFE densities at particle locations
for HACC-hydro particle data in GenericIO format

Requirements
------------

* GenericIO python interface: https://git.cels.anl.gov/hacc/genericio
* mpipartition: https://github.com/ArgonneCPAC/MPIPartition
* numpy

i.e.
```
pip install numpy mpipartition
pip install git+https://git.cels.anl.gov/hacc/genericio.git
```

Usage
-----

mpirun -n $PROCS python add_particle_dtfe_density.py <input> <RL> <OL> <output>
"""

import argparse
import numpy as np
import pydtfe
import pygio
from mpipartition import Partition, distribute, overload


def read_and_overload(filename, partition, rl, ol):
    data = pygio.read_genericio(filename)
    n_local = len(data["x"])

    data = distribute(partition, rl, data, ("x", "y", "z"))
    data["status"] = partition.rank * np.ones_like(data["x"], dtype=np.int32)
    data = overload(partition, rl, data, ol, ("x", "y", "z"))
    return data, n_local


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", type=str)
    parser.add_argument("RL", type=float)
    parser.add_argument("OL", type=float)
    parser.add_argument("output", type=str)
    parser.add_argument(
        "--species",
        type=str,
        default='bar',
        choices=['bar', 'all', 'dm', 'gas']
    )
    args = parser.parse_args()

    partition = Partition()
    if partition.rank == 0:
        print(f"Running with {partition.nranks} ranks")
        print(f"Decomposition: {partition.decomp}")
        print(
            f"Reading and Overloading Data (RL={args.RL}, OL={args.OL})...", flush=True
        )
    partition.comm.Barrier()

    data, n_local = read_and_overload(args.filename, partition, args.RL, args.OL)
    n_global = partition.comm.allreduce(n_local)

    if args.species == "bar":
        mask = np.bitwise_and(data["mask"], 4) > 0
    elif args.species == "gas":
        mask = (np.bitwise_and(data["mask"], 4) > 0) \
             & (np.bitwise_and(data["mask"], 8) == 0) \
             & (np.bitwise_and(data["mask"], 16) == 0) \
             & (np.bitwise_and(data["mask"], 128) == 0)
    elif args.species == "dm":
        mask = np.bitwise_and(data['mask'], 4) == 0
    elif args.species == "all":
        mask = np.ones_like(data["x"], dtype=np.bool_)
    else:
        raise NotImplementedError("unknown species")

    if partition.rank == 0:
        print("Done")
        print("Estimating particle volume and density...", flush=True)
    partition.comm.Barrier()

    pos = np.column_stack(
        [data["x"][mask], data["y"][mask], data["z"][mask]]
    )
    mass = data["mass"][mask][0]
    data["rho_dtfe"] = np.zeros_like(data["rho"])
    data["rho_dtfe"][mask] = pydtfe.compute_particle_density(pos, mass)

    partition.comm.Barrier()

    if partition.rank == 0:
        print("Done")
        print("Writing output...", flush=True)
    partition.comm.Barrier()

    # drop overloaded particles
    status = data.pop("status")
    mask_status = status == partition.rank
    for k in data.keys():
        data[k] = data[k][mask_status]

    # assert we have the correct number of global particles
    n_local_write = len(data["x"])
    n_global_write = partition.comm.allreduce(n_local_write)
    assert n_global == n_global_write

    pygio.write_genericio(
        args.output,
        data,
        [args.RL, args.RL, args.RL],
        [0.0, 0.0, 0.0],
        method=pygio.PyGenericIO.FileIO.FileIOMPICollective,
    )

    if partition.rank == 0:
        print("Done", flush=True)
    partition.comm.Barrier()
