# omniCADD
#### a CADD scoring system to asses variant deleteriousness in non-model organisms
----

This is where the updates for the CADD pipeline generalised to be able to acore more than one species will end up.
As it is still in its starting phase, it is privated and will be updated on the go.


## HPC Containers: Prefetch & Preflight
- Prefetch recommended images on a login node (avoids compute-node pulls):

```bash
# Augustus (+snpEff+Java) and optional snpEff image
bash scripts/prefetch_containers.sh
```

- Run a quick preflight to validate images and config before submitting jobs:

```bash
python3 scripts/preflight_check.py --profile config/test_overlay_augustus.yaml
# or for limited-data chr2 profile (no containers required):
python3 scripts/preflight_check.py --profile config/profile_limited_data.yaml
```

- Submit a small test on Dardel (includes preflight and Slurm profile):

```bash
sbatch scripts/dardel_submit_augustus_tests.slurm
```

- Local minimal chr2 run (preflight + conda):

```bash
bash scripts/run_hpc_chr2.sh
```

