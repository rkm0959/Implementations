# Heavy Ball method doesn't converge

In 2014, the paper by Lessard, Recht, and Packard showed (among many other incredible results) that Polyak's Heavy Ball method does not necessarily converge in the strongly-convex setting when we do not have `f` to be twice differentiable. 

Paper Link : https://arxiv.org/pdf/1408.3595.pdf

This repository aims to mimic their steps toward the proof. To be more precise, it does

- Using Weighted Off-By-One IQC to reproduce Figure 5 of the paper
- Implementing the author's counterexample and checking the method doesn't converge

Brief explanation of each file

- IQC.py : Implementation of Weighted Off-By-One IQC
- resultHeavy.PNG : Figure 5, reproduced
- counterexample.py : Implementation of the author's counterexample
- counterexample.PNG : The iterates of the heavy ball method which doesn't converge