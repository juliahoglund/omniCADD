Interactive Workflow Diagram
=============================

This page contains the full interactive workflow diagram for omniCADD.

.. raw:: html

   <style>
      body { margin: 0; padding: 0; }
      .document { max-width: 100% !important; }
   </style>
   <iframe src="_static/omniCADD_workflow_dag.html" 
           width="100%" 
           height="1150px" 
           style="border: none; display: block; margin: 0;">
   </iframe>

Features
--------

* **Hover** over any node to see configuration details
* **Toggle** Standard Tier, Limited Tier, and Conservation paths using the buttons
* **Color-coded lanes** show different pipeline stages:
  
  * Blue: Standard Tier (alignment-based workflow)
  * Purple: Limited Tier (alignment-free workflow) 
  * Green: Core Pipeline (shared steps)
  * Orange: Conservation Annotation (optional)
  * Yellow: Model Training & Output

* **Input nodes** (dashed borders) show required data files
* **Output nodes** (golden) highlight the final CADD scores

Direct Link
-----------

`Open diagram in full screen <_static/omniCADD_workflow_dag.html>`_ (opens in new window)

Pipeline Steps
--------------

For detailed descriptions of each step, see the :doc:`workflow` documentation.
