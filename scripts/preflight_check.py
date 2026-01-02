#!/usr/bin/env python3
import argparse
import os
import sys
try:
    import yaml
except Exception:
    yaml = None

REQUIRED_FILES = [
    # SNPEff DB inputs
    ("resources/snpEff/snpEff.config", True),
]

DERIVED_TEMPLATE = "results/derived_variants/singletons/chr{chr}.vcf"
GENOME_TEMPLATE = "resources/genome/Wild_Boar_chr{chr}.fa"


def check_exists(path: str, critical: bool = True):
    ok = os.path.exists(path)
    status = "OK" if ok else "MISSING"
    print(f"[{'CRIT' if critical else 'WARN'}] {status}: {path}")
    return ok or not critical


def main():
    ap = argparse.ArgumentParser(description="Preflight validation for omniCADD SNPEff run")
    ap.add_argument("--profile", required=True, help="YAML config file to read chromosomes and annotation paths from")
    args = ap.parse_args()

    if yaml is None:
        print("[CRIT] Missing dependency: pyyaml. Install with 'pip install pyyaml' or run with --profile pointing to JSON.")
        sys.exit(4)
    with open(args.profile, "r") as fh:
        cfg = yaml.safe_load(fh)

    # Check base required files
    ok = True
    for path, critical in REQUIRED_FILES:
        ok = check_exists(path, critical) and ok

    # Chromosomes
    chrs = cfg.get("chromosomes", {}).get("karyotype", [])
    if not chrs:
        print("[CRIT] No chromosomes defined in profile")
        sys.exit(2)

    # Gene annotation source
    gene_src = cfg.get("gene_annotation", {}).get("source", "gff")
    if gene_src == "gff":
        gff = cfg.get("gene_annotation", {}).get("gff", "")
        ok = check_exists(gff, True) and ok
    elif gene_src == "gtf":
        gtf = cfg.get("gene_annotation", {}).get("gtf", "")
        ok = check_exists(gtf, True) and ok
    else:
        print(f"[WARN] gene_annotation.source={gene_src} (no direct file check)")

    # Per-chromosome genome and derived VCFs
    for chr_ in chrs:
        genome = GENOME_TEMPLATE.format(chr=chr_)
        derived = DERIVED_TEMPLATE.format(chr=chr_)
        ok = check_exists(genome, True) and ok
        ok = check_exists(derived, False) and ok
        gz = derived + ".gz"
        tbi = gz + ".tbi"
        if os.path.exists(gz):
            ok = check_exists(tbi, True) and ok
        else:
            print(f"[WARN] Not gzipped/indexed: {derived}")

    print("Preflight status:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
