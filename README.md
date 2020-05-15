# Sandpiles

The aim of this repository is provide simulations of various **sandpile** models.
As of this writing the repository contains simulation of the **Abelian Sandpile** model introduced in <a href="#ref-BTW">[2]</a> (coded both in Julia and in C++ mainly for comparing the runtime speeds) and **boundary sandpile** model (coded in Julia) introduced in <a href="#ref-AS">[1]</a>.

## Abelian Sandpile

<a href = "https://en.wikipedia.org/wiki/Abelian_sandpile_model">Abelian Sandpile</a> was introduced in 1987 by <a href="#ref-BTW">Bak, Tang, and Wiesenfeld</a> as an example of a dynamical system exhibiting <a href = "https://en.wikipedia.org/wiki/Self-organized_criticality">self-organized criticality</a>.
The model, which is a cellular automaton, works according to the following simple rules. Start with `n` (indistinguishable) grains of  sands placed at the origin of the 2d integer lattice (i.e. the set of points of the 2-dimensional plane with integer coordinates). At any point of (discrete) time if a lattice site has at least 4 grains of sand, then it gives 1 grain to each of its 4 neighbours and loses 4 grains itself (this step is called `toppling`). This process of one-by-one topplings is being repeated until all lattice sites remain with at most 3 grains of sand. It can be proved that however large `n` might be, after finitely many toppling (the number depending on `n`) all lattice points will have at most 3 grains of sand. Moreover, it does not matter in which order the topplings are performed, the final configuration is always the same (hence the name **Abelian**). Unexpectedly, such a simple rule produces complex configurations. For example, the
final configuration of the Abelian sandpile for 1 000 000 grains of sand placed at the origin of 2d integer lattice looks as follows

<p align="center">
  <img src ="https://github.com/hayk314/Sandpiles/blob/master/C%2B%2B/AbelSand/Debug/Abel1000000.png" alt = "Abelian Sandpile">
</p>
<p align="center">
The final configuration of Abelian Sandpile on 2d grid for 1 000 000 particles at the origin. The coloring scheme (number of sand grains - color) is as follows: <b>0-black, 1-magenta, 2-red, 3-blue.</b>
</p>

### Running the code

#### C++ version

The complied code is the `/C++/AbelSand/Debug/AbelSand.exe`. Run it and follow its instructions, the user needs to fill the number of grains at the origin and hit the Enter button. The result is a **CSV** file showing the final distribution of the sandpile (`-1` there stands for the lattice site which was never visited by the process, in contrast to `0` meaning that vertex toppled at some point but eventually ended up with no grains in the final configuration). Along with the CSV file there is also a .ppm file (**ppm** stands for <a href ="https://en.wikipedia.org/wiki/Netpbm_format">portable pixmap format</a>) which can be converted into an image file using free online services such as <a href ="https://convertio.co/ppm-png/">here</a> or <a href ="https://www.freefileconvert.com/">here</a>.

Some benchmarks: on my `DELL notebook, Windows 10, Intel(R) Core(TM) i7-8550U CPU @ 1.80GHz (QuadCore)`, for 100 000 grains at the origin, 
the program runs in 8 seconds.


#### Julia version

The Julia program is located in `/Julia`. From Julia (version 1.1) `cd` to that folder and use
```Julia

A = include("Abel.jl")

Z_lat, Odometer = A.move(100000);
```
to run the program for 100 000 grains at the origin. This will produce the CSV files of the final configuration, the odometer function (the function on integer lattice counting for each lattice site how many times it has toppled until stabilization), and the colour image file of the final configuration (just like in the example image above).
On the same notebook as above, the program runs in 25 seconds.


As one can observe, the C++ version is faster than the Julia version, however not fast enough in general. While there are no heavy computations being performed, the process, however, does an excessive number of 2d array read-writes which eventually slow down the whole process. It would be very interesting to find out better approaches for stabilization than a straightforward simulation of the sandpile.

## Boundary sandpile

Boundary sandpile model was introduced by <a href="#ref-AS">Hayk Aleksanyan and Henrik Shahgholian</a> as an attempt to model a <em>qudrature surface</em> (a potential theoretic concept) by a particle dynamics. Contrary to the initial guess, however, the model turned out to produce a new and not yet fully understood phenomenon. 

To define the boundary sandpile model, assume we have mass `n>0` at the origin of the integer lattice (here `n` does not need to be an integer, it is a continuous mass). Fix also a threshold value equal to `sqrt(n)`, a boundary capacity of the model (this specific choice of the threshold is determined from scaling considerations, see <a href="#ref-AS">the original paper</a> for the details, in particular, for dimension `d >= 2` the value of the threshold is taken to be `n^{1/d}`). For each epoch of discrete time, define also a set of `visited sites` of the model, which are the set of all points of the integer lattice which have been visited by the spread of the process by the given time. For instance, at time `0` the only visited site would be the origin. At each point of time a given vertex of the lattice carrying a positive mass can topple if either it is in the interior of the set of visited sites, or it carries mass more than the `boundary threshold`. Toppling will distribute all the mass of the vertex equally among its neighours and will leave no mass for the vertex itself. Except for trivial cases, it takes an infinite number of topplings to stabilize the process. Nevertheless it is proved in <a href="#ref-AS">[A-S]</a> that once each vertex topples an infinite number of times the sandpile will always stabilize in the limit and the limiting configuration is independent of the order of the topplings (so the process is Abelian in this sense).

<p float ="center">
On the left is the final configuration of Boundary Sandpile on 2d grid for mass 1 000 at the origin, and the right image is the final configuration for mass 10 000. The colored pixels represent the boundary of the sandpile, where the entire mass is concentrated.
Lighter colors mean  more mass.
</p>

<p float="center">
 <span float = "left"> &nbsp; <img src="https://github.com/hayk314/Sandpiles/blob/master/Julia/BSand_Z_1000.png" width="300"  /> </span>
  <span>  &nbsp; &nbsp;  &nbsp;&nbsp;</span>
 <span> <img src="https://github.com/hayk314/Sandpiles/blob/master/Julia/BSand_Z_10000.png" width="500" />  </span>
</p>






## References

<li id="ref-AS">1. Hayk Aleksanyan, and Henrik Shahgholian  <a href = "https://arxiv.org/abs/1607.01525">Discrete Balayage and Boundary Sandpile</a>, Journal d'Analyse Mathematique (<i>to appear</i>) </li> 


<li id="ref-BTW">2. Per Bak, Chao Tang, and Kurt Wiesenfeld <a href = "https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.59.381">Self-organized criticality: An explanation of the 1/f noise</a>, Phys. Rev. A (3) 38, 1988</li>


<li id="ref-LPer">3. Lionel Levine, and Yuval Peres <a href = "https://arxiv.org/abs/1611.00411">Laplacian growth, sandpiles and scaling limits</a>, Bulletin (New Series) of the Amer. Math. Soc, 2017</li>


<li id="ref-LProp">4. Lionel Levine, and James Propp <a href ="https://www.ams.org/notices/201008/rtx100800976p.pdf">What is ... a Sandpile ?</a>, Notices Amer. Math. Soc, 2010</li>

