#!/usr/bin/env python3
import argparse
import os
import sys
try:
    import yaml
except Exception:
    yaml = None

DERIVED_TEMPLATE = "results/derived_variants/singletons/chr{chr}.vcf"


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

    ok = True

    # SNPEff config file from overlay/base config
    snpeff_conf = (
        cfg.get("annotation", {})
           .get("snpeff", {})
           .get("build", {})
           .get("config_file")
    )
    if snpeff_conf:
        ok = check_exists(snpeff_conf, True) and ok
    else:
        print("[CRIT] Missing annotation.snpeff.build.config_file in profile")
        ok = False

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
    genome_wildcard = (
        cfg.get("generate_variants", {})
           .get("reference_genome_wildcard")
        or "resources/genome/Wild_Boar_chr{chr}.fa"
    )
    for chr_ in chrs:
        try:
            genome = genome_wildcard.format(chr=chr_)
        except Exception:
            print(f"[CRIT] Invalid reference_genome_wildcard: {genome_wildcard}")
            sys.exit(3)
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
