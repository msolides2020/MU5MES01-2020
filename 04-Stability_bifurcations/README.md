#  Report 2
*MU5MES01 - 2019/20 - S.Brisard, D.Duhamel, C.Maurini, S. Neukirch*

The key concepts introduced in the second part of the class are

- Introduction to nonlinear elasticity
- FEniCS implementation of a non-linear elastic model
- Use and implementation of a Newton solver for the solution of nonlinear problems
- Solution of quasi-static rate-independent problems by time-stepping (progressive increment of a single loading parameter)
- Linearization of the nonlinear equilibrium equation within a variational formulation
- Linearized buckling analysis and solution of eigenvalue problems
- Stability analysis
- Discussion of the issues arising in the numerical solutions of nonlinear systems: possible non-convergence, convergence to an unstable solution, need for the introduction of imperfections, etc…


Your report should summarize and present synthetically your work on these items using the classical example of buckling of a column under its own weight as model problem. We give below some hints on how to write the report. Personalized analyses and comments are particularly welcome. Your are not obliged to follow the following format step by step. But you should include in your report the key concepts and results.


**Some suggestions:**

- Write concisely and effectively.
- Comment your results.
- The quality of the figures is important.
- Report only the minimal number of figures (of excellent quality) to effectively communicate your results.
- You can write in English or French.
- Use Latex for writing your report.
- In the written report you should correctly formulate each mathematical problem solved. You should not report all the details of the derivation of the formulation in the report. You will be asked about that during the oral examination.

**Important informations:**
  - Deadline: For the final version **Monday 4 January, 23h59**. 
  - **The maximal length of the report is 4 pages.**
  - To submit your report: 
      - An electronic version should be submitted to github. Proceed as follows to create the work and submission repository for your group:
        - **Only one of the two students** of your group will go to  https://classroom.github.com/g/N7o4GEQ4, accept the assignement and create a `new team`, naming the team as `NAMESTUDENT1-NAMESTUDENT2`.
        - **Once the first student has create the team**, the second student goes to https://classroom.github.com/g/N7o4GEQ4, accepts the assignment and asks to join the team with his name (do not create another team, there should one team for group).

      - In your "group" repository you should: 
          1. Create a directory called `CR2`
          2. Put your report in the pdf form named as `MES01-CR2-studentname1-studentname2.pdf` (file with a different naming scheme will not be accepted and evaluated. 
          3. Put all your files you used to obtain your results in `CR2/src` (namely the *.py and *.ipynb files)
  - We will evaluate the quality of the presentation (language, typesetting, and figures). Being able to effectively communicate your results is important for your future.
  - We ask you to be able to use git at least to push your data to the repository. This is the main reason why we ask to submit your report on the github platform. We will not accept submissions by mail.
# Nonlinear design of concrete bridge towers

The largest bridges in the world are [cable-stayed
bridges](https://en.wikipedia.org/wiki/Cable-stayed_bridge) and [suspension
bridges](https://en.wikipedia.org/wiki/Suspension_bridge), with spans that are
greater than 1000m. In both cases, the traffic loads and the deck self-weight
are transfered to the foundations through towers (or pylons) that are slender
structures, mainly subjected to compression. The figure below displays the
Millau cable-stayed bridge (credits: Nicolas Janberg, via
[Structurae](https://structurae.net/fr/medias/283009-viaduc-de-millau)).

![The Millau cable-stayed bridge (credits: Nicolas Janberg, via Structurae, https://structurae.net/fr/medias/283009-viaduc-de-millau).](283009-viaduc-de-millau.jpg)

Owing to their slenderness, these structures are subjected to so-called
*second-order effects*, which might lead to their ruin, just like compressing a
spaghetti might lead to its rupture. The goal of this project is to better
understand the role of the second-order effects.

We consider a simplified situation, where the tower is subjected to its own
weight only. In the remainder of this assignement, `x` denotes the direction of
gravity, `y` is the in-plane direction perpendicular to `x` and `z` is the
out-of-plane direction.

The pylon will be modelled as a hyperelastic, rectangular beam (length: `Lx`,
cross-section: `Ly × Lz`), clamped at one end and subjected to its own weight.

## Part 0. Introduction

The two questions below are based on strength considerations only. Stability
issues and second-order effects will be ignored. In other words, it will be
assumed that the deflection of the pylon under its own weight is negligible (the
pylon remains straight in its current configuration), see figure below.

**Question 0.1:** What is the maximum length `Lx` of a pylon made of concrete
(density: `ρ = 2400 kg/m³`; compressive strength: `fc = 50 MPa`)? Does that seem
realistic to you?

In reality, the geometry of the built structure is never perfect. In particular,
we assume here that the pylon, although perfectly straight, makes a small angle
`α` with the direction of gravity. (We use `α`=0.005 radians).

**Question 0.2:** Lack of verticality may introduce some tension in the
structure, and as concrete cannot sustain any tension, this will introduce
another limit in the length `Lx`.

Show that the maximum length `Lx` of a pylon is in fact smaller than the value
obtained before. To do so, apply elementary (linear) beam theory, observing that
the structure is statically determinate.

1. Compute the axial force and bending moment at any point of the pylon.
2. Compute the normal stress at any point of the pylon.
3. Show on a sketch of the pylon where the maximum compressive stress and maximum tensile
   stress occur.
4. Write two inequalities expressing the fact that the maximum compressive stress cannot exceed `fc`, and that the maximum tensile stress cannot be above 0.
5. Draw a graph in the plane `(Lx, Ly)` and show the two limiting curves, as well as the secure zone.
6. For a pylon with `Ly = 10m`, what is the maximum height `Lx`? What is the limiting cause in this case?
7. For a pylon with `Ly = 20m`, what is the maximum height `Lx`? What is the limiting cause in this case?
8. What should be the value of `Lz`?

![Second-order effects](second-order_effects.png)

## Second-order effects

The computations so far have not taken into account the deformation of the pylon
under the applied load. We will now take them into account and see that it lowers
dramatically the limiting values computed so far.

The second-order effects will be estimated numerically by means of a 2D (plane
stress) simulation using FEniCS. We consider the same beam as in the tutorial
`HyperelasticSolid.ipynb`. **Note that gravity refers to the abscissa in that
case!**

For the simulations, we will use the following numerical values

- beam length: `Lx = 1.0m`,
- Young modulus: `E = 1000Pa`,
- Poisson ratio: `ν = 0.3`.

We will test different values of the beam depth `Ly`, starting with `Ly = 0.1 * Lx`.

We define the following non-dimensional ratio

```
       ρ⋅g⋅Lx³
Γ = 12 ───────
        E⋅Ly²
```

When `α = 0` (pylon perfectly aligned with gravity), the straight configuration
becomes unstable when `Γ > 7.8`.

Above this critical value, the stable equilibrium configuration is deflected in
the `(x, y)` plane (see [Wikipedia](https://en.wikipedia.org/wiki/Self-buckling)).

Detecting this loss of stability requires the analysis of the eigenvalues of the
tangent stiffness matrix. A simpler approach is to introduce a geometric
imperfection that will trigger the instability and force the system to adopt its
buckled shape. We will use the angle `α` as an imperfection parameter to study
this buckling problem.

We suggest the following parametrization of the loading

```python
Gamma = 7.84
rho_g = Gamma * Y * Ly**2 / 12 / Lx**3
b0 = dolfin.Expression(("-t * rho_g * cos(alpha)", "t * rho_g * sin(alpha)"),
                       t=0, alpha=alpha, rho_g=rho_g, degree=0)
```

**Question 0.3:** what is the meaning of `t` ?

## Part I. Nonlinear quasi-static response analysis

Here we use

- beam length: `Lx = 1.0m`,
- Young modulus: `E = 1000Pa`,
- Poisson ratio: `ν = 0.3`.

The goal of this part is to reproduce a bifurcation diagram similar to the one
in the figure below where we report the transverse-end displacement `u_L / Lx`
as a function of the loading-parameter `t`

![schema](./equilibrium_curve.png)

**Question I.1:** Set `Ly = 0.1 * Lx`, `α = 0.1` and plot the normalized deflection
`u_L / Lx` at the free end as a function of `t`, for `0 ≤ t ≤ 2`. Interpret the
results.

**Question I.2:** Discuss the influence of the imperfection `α` on the diagram:
for example how the bifurcation curve is affected by `α`.

**Question I.3:** Explicitly present the algorithm you used (mesh, time-stepping and
Newton solver), and justify the value of parameters you used to perform the
simulation (load-step, minimum and maximum loads, numerical tolerances, initial
guesses, imperfection…).
The stability limit of the fundamental branch coincides with the first buckling load. Why? Are you able to prove it?

## Part II. Stability analysis

By studying the sign of the smallest eigenvalue of the hessian matrix, study the stability of the three branches in the bifurcation diagram given in part I. Produce a figure where each branch is coloured in red if unstable and in green if stable. The stability limit of the fundamental branch coincides with the first buckling load. Why? Are you able to prove it?

For the stability analysis, you can use the notebook `NonlinearStabilityCheck.ipynb` as an example 
Further documentation:
  - Some details on the linearization of the equilibrium equations and the derivatives can be found in Section 3.5 of [2]
  - You can look also to [3] for the buckling analysis.


**Advanced (to do only if you have finished the rest)**
-  Estimate the buckling load for the vertical configuration by performing a linearized buckling analysis by solving and eigenvalue problem of the type`(K-\lambda G) U =0` where `K` and `G` are the matrices obtained by assembling the bilinear forms associated to the elastic and geometric stiffness, respectively (see the pdf of the presention 11/12/2020). 

*Further documentation*

- Some details on the linearization of the equilibrium equations and the
  derivatives can be found in Section 3.5 of [2]
- You can look also to [3] for the buckling analysis.

## Part III. Going further

We now use parameters value for concrete `g = 9.81 (SI)`, `E = 35GPa`, `α =
0.005rad`, density `ρ = 2400kg/m³`, and compressive strength `fc = 50MPa`.

**Question III.1:** For a pylon with base `Ly = 10m`, what is the height `Lx`
for which we reach `Γ = 7.8`?

**Question III.2:** For a pylon with base `Ly = 10m`, what is the maximum height
`Lx`?  Is this limition due to tension or compression? Compare to the value
found in Part I. Comment.

**Question III.3:** For a pylon with base `Ly = 20m`, what is the height `Lx`
for which we reach `Γ = 7.8`?

**Question III.4:** For a pylon with base `Ly = 20m`, what is the height `Lx`
for which the limitation is due to tension? Due to compression? Compare to the
value found in Part I. Comment.

## Part 999. Only if you dont take vacations

We note `P_i` the different bifurcation points.

Here we use

- beam length: `Lx = 1.0m`, `Ly = 0.1m`
- Young modulus: `E = 1000Pa`,
- Poisson ratio: `ν = 0.3`.

**Question 999.1:** Find the second bifurcation point `P_2`, for example in the
following way. Go through `P_1` by setting `α = 0`, then use a small value of
`α` to hit the second bifurcation point `P_2`.

**Question 999.2** Is it true that for `Ly > 0.3` only `P_1` and `P_2` remain?

**Question 999.3** Is it ture that for `Ly > 0.9` there is no longer any bifurcation point?

## References
You find a full analytical solution of the buckling and post-buckling problem obtained when using a one-dimensional beam model in [1]

[1] J.-J. Marigo, Mécanique des Milieux Continus I, Notes de cours MEC 430,  Ecole Polytechnique, https://cel.archives-ouvertes.fr/cel-01023392

[2] P. Wriggers, Nonlinear finite element methods, Springer 2008

[3] F.Voldoire, Y.Bamberger, Mécanique des structures, Presse de l'Ecole Nationale des Ponts et Chaussées, 2008, chapitre 9, pag 521.

<!-- Local Variables: -->
<!-- fill-column: 80 -->
<!-- End: -->
