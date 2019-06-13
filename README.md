# Sandpiles

The aim of this repository is provide simulations of various **sandpile** models.
As of this writing the repository contains simulation of the **Abelian Sandpile** model introduced in <a href="#ref-BTW">[2]</a> (coded both in Julia and in C++ mainly for comparing the runtime speeds) and **boundary sandpile** model (coded in Julia) introduced in <a href="#ref-AS">[1]</a>.

## Abelian Sandpile

<a href = "https://en.wikipedia.org/wiki/Abelian_sandpile_model">Abelian Sandpile</a> was introduced in 1987 by <a href="#ref-BTW">Bak, Tang, and Wiesenfeld</a> as an example of a dynamical system exhibiting <a href = "https://en.wikipedia.org/wiki/Self-organized_criticality">self-organized criticality</a>.
The model, which is a cellular automaton, works according to the following simple rules. Start with `n` (indistinguishable) grains of  sands placed at the origin of the 2d integer lattice (i.e. the set of points of the 2-dimensional plane with integer coordinates). At any point of (discrete) time if a lattice site has at least 4 grains of sand, then give 1 grain to each of its 4 neighbours (this step is called `toppling`). Repeat this process one-by-one, until all sites of the lattice have at most 3 grains of sand. It can be proved that however large `n` might be, after finitely many toppling (the number depending on `n`) all lattice points will have at most 3 grains of sand. Moreover, it does not matter in which order the topplings are performed, the final configuration is always the same (hence the name **Abelian**). Unexpectedly, such a simple rule produces complex configurations. For example, the
final configuration of the Abelian sandpile for 1 000 000 grains of sand placed at the origin of 2d integer lattice looks as follows

<p align="center">
  <img src ="https://github.com/hayk314/Sandpiles/blob/master/C%2B%2B/AbelSand/Debug/Abel1000000.png" alt = "Abelian Sandpile">
</p>
<p align="center">
The final configuration of Abelian Sandpile on 2d grid for 1 000 000 particles at the origin. The coloring scheme (number of sand grains - color) is as follows: <b>0-black, 1-magenta, 2-red, 3-blue.</b>
</p>

## References

<li id="ref-AS">1. Hayk Aleksanyan, and Henrik Shahgholian  <a href = "https://arxiv.org/abs/1607.01525">Discrete Balayage and Boundary Sandpile</a>, Journal d'Analyse Mathematique (<i>to appear</i>) </li> 


<li id="ref-BTW">2. Per Bak, Chao Tang, and Kurt Wiesenfeld <a href = "https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.59.381">Self-organized criticality: An explanation of the 1/f noise</a>, Phys. Rev. A (3) 38, 1988</li>


<li id="ref-LPer">3. Lionel Levine, and Yuval Peres <a href = "https://arxiv.org/abs/1611.00411">Laplacian growth, sandpiles and scaling limits</a>, Bulletin (New Series) of the Amer. Math. Soc, 2017</li>


<li id="ref-LProp">4. Lionel Levine, and James Propp <a href ="https://www.ams.org/notices/201008/rtx100800976p.pdf">What is ... a Sandpile ?</a>, Notices Amer. Math. Soc, 2010</li>

