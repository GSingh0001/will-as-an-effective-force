# Observational Coupling Tests and Will as an Effective Force

Two companion preprints by Ganga Singh (Independent Researcher, ORCID [0009-0004-2191-356X](https://orcid.org/0009-0004-2191-356X)), with the frozen test harness behind the first paper. Both papers are published on Zenodo under CC-BY 4.0 and cross-linked there. This repository holds the same texts in markdown and PDF, plus the code, so the work can be read, indexed, and rerun without going through a PDF download.

## Paper 1: Born, Energy and No-Signalling Tests of Observational Coupling Models

A reproducible falsification study. Before a candidate coupling term added to quantum dynamics is interpreted, it should pass three cheap numerical tests: (K1) recovery of Born rule outcome frequencies from individual conditional trajectories, (K2) an energy audit consistent with the role claimed for the coupling, and (K3) a finite ensemble-consistency test of the kind motivated by Gisin's no-signalling argument. Standard continuous measurement passes all three. A deterministic Hermitian coupling fails K1 exactly. A deterministic state-dependent pull fails K3. Four further frozen-run candidates fail K1. A stochastic repair of the failed coupling is shown, by an explicit parameter substitution stated in the paper's appendix, to be identical to standard continuous measurement theory (Wiseman and Milburn; Belavkin), so its novelty question is closed.

- Zenodo: https://doi.org/10.5281/zenodo.24401299
- Text: [papers/paper1-kill-tests.md](papers/paper1-kill-tests.md) · PDF: [papers/Born-Energy-and-No-Signalling-Tests.pdf](papers/Born-Energy-and-No-Signalling-Tests.pdf)

## Paper 2: Will as an Effective Force

Units, channel selection coupling, and a thermodynamic floor. Will is defined operationally as a shift in the log odds of an outcome being selected, measured in nats against an unforced baseline. The paper tests the unit against metrological criteria, excludes outcome bias and noise steering as coupling points by the kill tests of Paper 1, identifies measurement-strength modulation as known continuous measurement physics, and keeps one surviving coupling: selection among channels and measurement bases, which leaves the statistics inside a selected channel unchanged. The thermodynamic section states a corrected accounting: under a stated Boltzmann implementation a log-odds shift is an energy-gap change of W times k_B T, and preparing a biased distribution is bounded below by k_B T times its relative entropy to baseline. A stronger earlier claim, a universal charge of k_B T per nat, is withdrawn in the text. The paper also records that its choice-domain dilution law is equivalent to divisive normalisation and loses out of sample to nested logit models, and it states its falsifier in the abstract.

- Zenodo: https://doi.org/10.5281/zenodo.24401500
- Text: [papers/paper2-effective-theory.md](papers/paper2-effective-theory.md) · PDF: [papers/Will-as-an-Effective-Force.pdf](papers/Will-as-an-Effective-Force.pdf)

## The harness

`killtests/` is the frozen reference implementation behind Paper 1: a single-qubit trajectory simulator and test runner (Euler-Maruyama, dt = 0.005, numpy PCG64, frozen seeds 20261118 to 20261122), with the harness report and the machine-readable results of the frozen runs.

- `run_killtests.py` and `trajectories.py`: the implementation.
- `KILL-TEST-REPORT.md`: protocols, acceptance criteria, and verdicts.
- `results.json`: the frozen run outputs the paper reports.

Every verdict in Paper 1 was reproduced in an independent regrade from byte-identical copies of this code on a different machine; Section 3.3 of the paper also gives a protocol for reproducing the verdicts from the paper alone.

## What this work does not claim

No candidate tested produced a discrepancy with standard quantum mechanics: every failure was a failure to reproduce it, and every pass was a reproduction of it. Paper 2's two surviving predictions (a capacity ceiling of 2.15 nats under incentive scaling, and veto dispatch collapse under channel load with an invariant point of no return) are stated as preregistered protocols awaiting human data; neither has been tested on people. The papers say this in their own text, and this README does not upgrade them.

## Licence

The code in `killtests/` is released under the MIT License (see [LICENSE](LICENSE)). The paper texts are published under CC-BY 4.0 on Zenodo and are reproduced here under the same licence.

## Citation

Singh, G. (2026). Born, Energy and No-Signalling Tests of Observational Coupling Models: A Reproducible Falsification Study. Zenodo. https://doi.org/10.5281/zenodo.24401299

Singh, G. (2026). Will as an Effective Force: Units, Channel Selection Coupling, and a Thermodynamic Floor. Zenodo. https://doi.org/10.5281/zenodo.24401500
