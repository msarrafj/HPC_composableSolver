<center> <h1> Codes for composable block solver methodologies for the four-field double porosity/permeability model </h1> </center>

Codes for 
><i>Mohammad S. Joshaghani, Justin Chang, Kalyana B. Nakshatrala, and Matthew G. Knepley</i>, "On composable block solvers and performance spectrum analysis for double porosity/permeability model" [Journal of Computational Physics](https://www.sciencedirect.com/science/article/pii/S0021999119301378), 386: 428-466, 2019. [Available in [arXiv](https://arxiv.org/abs/1808.08328)]
>
><details><summary>[Abstract]</summary>
><p> 
>The objective of this paper is twofold. First, we propose two composable block solver methodologies to solve the discrete systems that arise from finite element discretizations of the double porosity/permeability (DPP) model. The DPP model, which is a four-field mathematical model, describes the flow of a single-phase incompressible fluid in a porous medium with two distinct pore-networks and with a possibility of mass transfer between them. Using the composable solvers feature available in PETSc and the finite element libraries available under the Firedrake Project, we illustrate two different ways by which one can effectively precondition these large systems of equations. Second, we employ the recently developed performance model called the Time-Accuracy-Size (TAS) spectrum to demonstrate that the proposed composable block solvers are scalable in both the parallel and algorithmic sense. Moreover, we utilize this spectrum analysis to compare the performance of three different finite element discretizations (classical mixed formulation with H(div) elements, stabilized continuous Galerkin mixed formulation, and stabilized discontinuous Galerkin mixed formulation) for the DPP model. Our performance spectrum analysis demonstrates that the composable block solvers are fine choices for any of these three finite element discretizations. Sample computer codes are provided to illustrate how one can easily implement the proposed block solver methodologies through PETSc command line options. 
></p>
></details>

This repository provides computer codes for two recently proposed composable block solver methodologies to solve the discrete systems that arise from double porosity/permeability (DPP) model. 
Time-Accuracy-Size (TAS) spectrum performance model is employed to (i) demonstrate that solvers are scalable in both the parallel and algorithmic sense (ii) compare the performance of three different finite element discretizations (classical mixed formulation with H(div) elements, stabilized continuous Galerkin mixed formulation, and stabilized discontinuous Galerkin mixed formulation) for the DPP model. More details are discussed in the paper.

