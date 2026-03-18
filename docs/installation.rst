Installation
============

This guide covers installing and setting up omniCADD.

Prerequisites
-------------

System Requirements
~~~~~~~~~~~~~~~~~~~

* **Operating System**: Linux or macOS (Windows via WSL2)
* **RAM**: 16+ GB (64+ GB recommended for whole genomes)
* **Storage**: 500+ GB free space
* **CPU**: Multi-core processor (8+ cores recommended)

Software Requirements
~~~~~~~~~~~~~~~~~~~~~

Required Software
^^^^^^^^^^^^^^^^^

* **Snakemake** ≥ 7.0
* **Conda** or **Mamba** (mamba recommended for faster installation)
* **Git** for cloning the repository

Optional Software
^^^^^^^^^^^^^^^^^

* **Docker** or **Singularity/Apptainer** (for containerized execution)
* **SLURM** job scheduler (for HPC cluster execution)

Installation Steps
------------------

1. Clone the Repository
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone https://github.com/yourusername/omniCADD.git
   cd omniCADD

2. Install Conda/Mamba
~~~~~~~~~~~~~~~~~~~~~~

If you don't have conda installed:

**Miniconda** (recommended):

.. code-block:: bash

   # Linux
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   
   # macOS
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
   bash Miniconda3-latest-MacOSX-x86_64.sh

**Install Mamba** (optional but recommended):

.. code-block:: bash

   conda install -n base -c conda-forge mamba

3. Configure Conda Channels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Set up bioconda and conda-forge channels:

.. code-block:: bash

   conda config --add channels defaults
   conda config --add channels bioconda
   conda config --add channels conda-forge
   conda config --set channel_priority strict

4. Install Snakemake
~~~~~~~~~~~~~~~~~~~~

Create a dedicated environment for Snakemake:

.. code-block:: bash

   # With mamba (faster)
   mamba create -n snakemake -c conda-forge -c bioconda snakemake
   
   # Or with conda
   conda create -n snakemake -c conda-forge -c bioconda snakemake
   
   # Activate environment
   conda activate snakemake

5. Verify Installation
~~~~~~~~~~~~~~~~~~~~~~~

Check that Snakemake is installed correctly:

.. code-block:: bash

   snakemake --version
   # Should output: 7.x.x or higher

Test the workflow:

.. code-block:: bash

   # Dry run to check configuration
   snakemake -n

Optional: Container Setup
--------------------------

Docker Installation
~~~~~~~~~~~~~~~~~~~

**Linux**:

.. code-block:: bash

   # Install Docker Engine
   curl -fsSL https://get.docker.com -o get-docker.sh
   sudo sh get-docker.sh
   
   # Add user to docker group
   sudo usermod -aG docker $USER
   # Log out and back in for this to take effect

**macOS**:

Download and install Docker Desktop from https://www.docker.com/products/docker-desktop

Singularity/Apptainer (HPC Clusters)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most HPC systems have Singularity/Apptainer pre-installed. Check with:

.. code-block:: bash

   singularity --version
   # or
   apptainer --version

If not available, contact your system administrator.

Environment Setup
-----------------

The pipeline uses multiple conda environments for different steps. These are automatically 
created by Snakemake when you run the workflow with ``--use-conda``.

**Important**: The first run will take longer (30-90 minutes) as conda environments are created.

Pre-create Environments (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To pre-create all environments before running the pipeline:

.. code-block:: bash

   # Create all environments
   snakemake --use-conda --conda-create-envs-only --cores 1

Verify Environment Creation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check that environments were created:

.. code-block:: bash

   ls .snakemake/conda/
   # Should show multiple environment directories

Next Steps
----------

After successful installation:

1. ✅ Read :doc:`input_data` to prepare your genomic resources
2. ✅ Configure the pipeline in :doc:`configuration`
3. ✅ Run the workflow following the user guide

Additional Resources
--------------------

* **Snakemake Documentation**: https://snakemake.readthedocs.io/
* **Conda User Guide**: https://docs.conda.io/projects/conda/en/latest/user-guide/
* **Bioconda**: https://bioconda.github.io/
