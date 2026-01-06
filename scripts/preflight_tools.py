#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os
import shutil
import shlex
try:
    import yaml
except Exception:
    yaml = None

def run(cmd: list[str]) -> tuple[int, str]:
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
        return 0, out
    except subprocess.CalledProcessError as e:
        return e.returncode, e.output

BASE_DIR = os.getcwd()

def exec_in(img: str, args: list[str]) -> tuple[int, str]:
    # Resolve image path relative to the profile directory so relative paths work
    if not img.startswith("docker://") and not os.path.isabs(img):
        img = os.path.abspath(os.path.join(BASE_DIR, img))
    # Always execute in a login bash, cd into base dir, load modules, then run via apptainer or singularity.
    quoted_args = " ".join(shlex.quote(a) for a in args)
    cmd = (
        f"cd {shlex.quote(BASE_DIR)}; "
        "ml PDC >/dev/null 2>&1 || true; "
        "ml apptainer >/dev/null 2>&1 || true; "
        "if command -v apptainer >/dev/null 2>&1; then "
        f"apptainer exec {shlex.quote(img)} {quoted_args}; "
        "elif command -v singularity >/dev/null 2>&1; then "
        f"singularity exec {shlex.quote(img)} {quoted_args}; "
        "else echo 'Container runtime not found. Load modules: ml PDC; ml apptainer' >&2; exit 127; fi"
    )
    return run(["bash", "-lc", cmd])

def check_augustus(img: str) -> bool:
    ok = True
    code, out = exec_in(img, ["augustus", "--version"])
    print("[AUG] version rc=", code)
    if code != 0:
        print(out)
        return False
    # detect config path
    for p in [
        "/opt/conda/envs/augustus/share/augustus/config",
        "/opt/conda/share/augustus/config",
        "/usr/share/augustus/config",
        "/usr/local/share/augustus/config",
    ]:
        code, _ = exec_in(img, ["bash", "-lc", f"[ -d '{p}' ]"])
        if code == 0:
            print("[AUG] config:", p)
            break
    else:
        print("[AUG] WARN: could not find config dir")
    return ok

def check_snpeff(img: str) -> bool:
    ok = True
    code, out = exec_in(img, ["bash", "-lc", "which snpEff && snpEff -version"]) 
    print("[SNP] version rc=", code)
    if code != 0:
        print(out)
        ok = False
    code, out = exec_in(img, ["bash", "-lc", "java -version"]) 
    print("[SNP] java rc=", code)
    if code != 0:
        print(out)
        ok = False
    return ok

def main():
    ap = argparse.ArgumentParser(description="Preflight check for containers (Augustus, SNPEff)")
    ap.add_argument("--profile", required=True, help="Config YAML with containers and snpeff settings")
    args = ap.parse_args()
    if yaml is None:
        print("ERROR: pyyaml not installed", file=sys.stderr)
        sys.exit(2)
    profile_path = os.path.abspath(args.profile)
    global BASE_DIR
    BASE_DIR = os.path.dirname(profile_path)
    with open(profile_path, "r") as fh:
        cfg = yaml.safe_load(fh)
    cont = cfg.get("containers", {})
    aug = cont.get("augustus")
    sne = cont.get("snpeff")
    ok = True
    if aug:
        print("Checking Augustus container:", aug)
        ok = check_augustus(aug) and ok
    if sne:
        print("Checking SNPEff container:", sne)
        ok = check_snpeff(sne) and ok
    # SNPEff config file check
    snpeff_conf = (
        cfg.get("annotation", {})
           .get("snpeff", {})
           .get("build", {})
           .get("config_file")
    )
    if snpeff_conf:
        exists = os.path.exists(snpeff_conf)
        print(f"[SNP] config exists: {exists} -> {snpeff_conf}")
        ok = exists and ok
    if not ok:
        print("Preflight status: FAIL", file=sys.stderr)
        print("Hint: load modules first with 'ml PDC; ml apptainer' if missing.", file=sys.stderr)
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
