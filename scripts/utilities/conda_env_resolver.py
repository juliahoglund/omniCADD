#!/usr/bin/env python3
"""
Smart conda environment path resolver for omniCADD pipeline.
Automatically switches between generic and cluster-specific environments.
"""

import os
import sys
import socket

def detect_cluster():
    """Detect which cluster/environment we're running on."""
    hostname = socket.gethostname().lower()
    
    # Check for UPPMAX clusters (Rackham, Snowy, Bianca)
    if any(cluster in hostname for cluster in ['rackham', 'snowy', 'bianca']):
        return 'uppmax'
    
    # Check for PDC Dardel cluster
    if 'dardel' in hostname or 'pdc.kth.se' in hostname:
        return 'dardel'
    
    # Check for other common HPC indicators
    if any(indicator in hostname for indicator in ['login', 'compute', 'node']):
        return 'hpc'
        
    # Check environment variables that indicate HPC systems
    if any(var in os.environ for var in ['SLURM_JOB_ID', 'PBS_JOBID', 'LSB_JOBID']):
        return 'hpc'
    
    # Check for specific UPPMAX environment variables
    if 'UPPMAX_CLUSTER' in os.environ or 'SNIC_RESOURCE' in os.environ:
        return 'uppmax'
        
    # Check for PDC environment variables
    if 'PDC_CLUSTER' in os.environ or any(var.startswith('PDC_') for var in os.environ):
        return 'dardel'
        
    return 'local'

def get_conda_env_path(env_name):
    """Get the appropriate conda environment path based on the detected cluster."""
    cluster = detect_cluster()
    
    if cluster == 'uppmax':
        # Use Rackham-specific environments
        env_path = f"config/slurm/rackham/{env_name}.yml"
        if os.path.exists(env_path):
            return env_path
        else:
            print(f"Warning: Rackham-specific env not found: {env_path}", file=sys.stderr)
            print(f"Falling back to generic environment", file=sys.stderr)
    
    elif cluster == 'dardel':
        # Use Dardel-specific environments
        env_path = f"config/slurm/dardel/{env_name}.yml"
        if os.path.exists(env_path):
            return env_path
        else:
            print(f"Warning: Dardel-specific env not found: {env_path}", file=sys.stderr)
            print(f"Falling back to generic environment", file=sys.stderr)
    
    # Default to generic environments
    return f"workflow/envs/{env_name}.yml"

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python conda_env_resolver.py <env_name>", file=sys.stderr)
        sys.exit(1)
        
    env_name = sys.argv[1]
    print(get_conda_env_path(env_name))
