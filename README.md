# MSMPC: Multi-Subset Multi-Particle Correlators

Simple code to automatically compute n-particle correlators between n or less arbitrary subsets (n=2,…∞), without the need of lengthy analytic formulae. Based on the recursive algorithm from the "Generic Framework" (arXiv:1312.3572).

Written in C++, from AliRoot/AliPhysics classes.

Running example based on the toy MC "FlowAnalysisOnTheFly" (credits: Ante Bilandzic), execute with:

$ root runFlowAnalysisOnTheFly.C

Contact: Jacopo Margutti (jacopo.margutti@cern.ch)

***********************************************************************************************************************************

Notes:

1) Subsets are currently defined only according to particles' charge and pseudorapidity, can be trivially extended to any particle feature.
2) The recursive algorithm to compute correlators is encoded in the method AliFlowAnalysisSubsetsGF::ucN()



