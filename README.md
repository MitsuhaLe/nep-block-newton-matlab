# BlockNewtonNEP

A MATLAB/Julia reproduction project for the block Newton method proposed in

> D. Kressner,  
> *A Block Newton Method for Nonlinear Eigenvalue Problems*,  
> Numerische Mathematik, 114, 355–372 (2009).
> https://doi.org/10.1007/s00211-009-0259-x
> 
This repository reproduces the numerical experiments from the paper and
includes additional implementations of the `gun` problem in both Julia and MATLAB.

One motivation of this project is to provide a MATLAB implementation
corresponding to the Julia implementation in NEP-PACK,
making the block Newton method easier to study, debug, and verify
for MATLAB users.

A **detailed reproduction report** is available in:

```text
reproduction_report.pdf
```

---

## Features

- Reproduction of the numerical experiments in Kressner's paper
- MATLAB implementation of the Block Newton method
- Julia implementation for comparison and verification
- Additional `gun` problem example
- Cross-validation between Julia and MATLAB implementations

---

## Project Structure

```text
BlockNewtonNEP/
├─ ex13_14_15/                # MATLAB scripts for Examples 13--15
│  ├─ example13.m
│  ├─ example14.m
│  ├─ example15.m
│  ├─ armijoRule.m
│  ├─ initialPair.m
│  ├─ nlevp_newtonstep.m
│  └─ resT.m
│
├─ testgun_julia/             # Julia implementation of the gun problem
│  └─ testgun.jl
│
├─ testgun_matlab/            # MATLAB implementation of the gun problem
│  ├─ testgun.m
│  ├─ gun.m
│  ├─ armijoRule.m
│  ├─ nlevp_newtonstep.m
│  ├─ resT.m
│  └─ gun_curve.eps
│
├─ others/                    # Simple LU factorization timing comparison
│  ├─ lutime.jl
│  └─ lutime.m
│
└─ method_blocknewton.jl      # Original NEP-PACK implementation
```

---

## Results

example13: The curves have similar shapes.

<table>
	<tr>
		<th>Original curve</th>
		<th>Reproduced curve</th>
	</tr>
	<tr>
		<td><img src="/others/ex13origin.png" alt="original curve" /></td>
		<td><img src="/others/ex13res.png" alt="reproduced curve" /></td>
	</tr>
</table>

---

## Quick Start

### Clone the repository

```bash
cd BlockNewtonNEP
```

---

## MATLAB

Run the numerical examples:

```matlab
example13
example14
example15
```

Run the `gun` problem example:

```matlab
cd testgun_matlab
testgun
```

---

## Julia

Run the Julia version of the `gun` problem:

```julia
include("testgun_julia/testgun.jl")
```

---

## Numerical Experiments

The repository reproduces the numerical experiments presented in:

- Example 13
- Example 14
- Example 15

from Kressner's paper, together with an additional `gun` problem example.

The Julia and MATLAB implementations can be compared directly
to verify correctness and consistency.

---

## Dependencies

- Julia: NEP-PACK
- MATLAB: NLEVP 

---

## References

[1] D. Kressner,  
*A Block Newton Method for Nonlinear Eigenvalue Problems*,  
Numerische Mathematik, 114, 355–372 (2009).

[2] NEP-PACK  
https://github.com/NEP-Pack/NEP-Pack

[3] NLEVP Collection  
https://github.com/ftisseur/nlevp

[4] Julia Language  
https://github.com/JuliaLang/julia

---

## License

MIT License.

---

## Acknowledgments

This project was inspired by the implementation in NEP-PACK
and the nonlinear eigenvalue problem collection NLEVP.

Contributions, issues, and suggestions are welcome.