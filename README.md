**Numerical study and detection of Resonances with Python**

Resonances are essentially particles that can decay into other particles. 
They can be defined as poles (singularities) in a complex plane defined by the energy and the inverse of the lifetime of the particle. 
In Hadron Physics these poles can be simulated numerically using Effective Theories. 

This project consists on a set of programs used to optimize/automatize all things involved in the description of these resonances/poles.

This project consists on the following python programs:
- **detectmaxima.py**: Main program.
This is a specialized program that detects all the poles in a given region, plotting also the thresholds. It can detect poles in several Riemann sheets. Later versions calculate the couplings and the contribution of each channels to the resonance wave function (compositeness).
- **matrix.py**:
This is an input program where we add the interaction matrices, masses and other parameters needed in functions.py.
- **functions.py**:
In this program we add the interaction potential and the loop funcion(s) than will then be used to calculate the T-matrix.
- **polescalc.py**:
This program takes the inputs of matrix.py and functions.py and calculates the T-matrix, as well as other auxiliary python functions needed to plot the amplitudes in the energy complex plane.
- **poles.py**:
Here we have the pole class, which takes as inputs the mass and width of the pole (an initial guess can be obtained from the detectmaxima.py program) as well as any other necessary parameters. It has methods that can find the pole position, calculate the couplings, follow the poles while varying any of the parameters, compare the couplings with any other pole calculating a distance (for example the euclidian distance - this can be used to compare poles with different parameters that may actually be the same pole), plot the region in the complex plane around the pole.

This was one of the main programs I used to complete my PhD. 
Thesis: https://inspirehep.net/literature/1827249
