# Boussinesq thermal convection in a spherical shell

The model equations are

<img src="https://latex.codecogs.com/svg.latex?\inline&space;\begin{align*}&space;\left(\partial_t&space;-&space;\Delta\right)\mathbf{u}&space;&&space;=&space;\mathbf{u}&space;\times&space;\left(&space;\nabla&space;\times&space;\mathbf{u}\right)&space;&plus;&space;\frac{Ra}{Pr}&space;\Theta&space;\mathbf{r}&space;-&space;\nabla\Pi&space;\\[0.6cm]&space;\left(\partial_t&space;-&space;\frac{1}{Pr}\Delta\right)\Theta&space;&&space;=&space;\frac{S}{Pr}&space;-&space;\mathbf{u}\cdot\nabla\Theta\\[0.6cm]&space;\nabla\cdot\mathbf{u}&space;&&space;=&space;0&space;\end{align*}" title="\begin{align*} \left(\partial_t - \Delta\right)\mathbf{u} & = \mathbf{u} \times \left( \nabla \times \mathbf{u}\right) + \frac{Ra}{Pr} \Theta \mathbf{r} - \nabla\Pi \\[0.6cm] \left(\partial_t - \frac{1}{Pr}\Delta\right)\Theta & = \frac{S}{Pr} - \mathbf{u}\cdot\nabla\Theta\\[0.6cm] \nabla\cdot\mathbf{u} & = 0 \end{align*}" />
![Model equations](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cbegin%7Balign*%7D%20%5Cleft%28%5Cpartial_t%20-%20%5CDelta%5Cright%29%5Cmathbf%7Bu%7D%20%26%20%3D%20%5Cmathbf%7Bu%7D%20%5Ctimes%20%5Cleft%28%20%5Cnabla%20%5Ctimes%20%5Cmathbf%7Bu%7D%5Cright%29%20&plus;%20%5Cfrac%7BRa%7D%7BPr%7D%20%5CTheta%20%5Cmathbf%7Br%7D%20-%20%5Cnabla%5CPi%20%5C%5C%5B0.6cm%5D%20%5Cleft%28%5Cpartial_t%20-%20%5Cfrac%7B1%7D%7BPr%7D%5CDelta%5Cright%29%5CTheta%20%26%20%3D%20%5Cfrac%7BS%7D%7BPr%7D%20-%20%5Cmathbf%7Bu%7D%5Ccdot%5Cnabla%5CTheta%5C%5C%5B0.6cm%5D%20%5Cnabla%5Ccdot%5Cmathbf%7Bu%7D%20%26%20%3D%200%20%5Cend%7Balign*%7D)

with the parameters defined as

<img src="https://latex.codecogs.com/svg.latex?\inline&space;\begin{align*}&space;Pr&space;&&space;=&space;\frac{\nu}{\kappa}\\&space;Ra&space;&&space;=&space;\frac{g&space;\alpha&space;\beta&space;r_o^4}{\nu\kappa}&space;\end{align*}" title="\begin{align*} Pr & = \frac{\nu}{\kappa}\\ Ra & = \frac{g \alpha \beta r_o^4}{\nu\kappa} \end{align*}" />
![Nondimensional parameters](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cbegin%7Balign*%7D%20Pr%20%26%20%3D%20%5Cfrac%7B%5Cnu%7D%7B%5Ckappa%7D%5C%5C%20Ra%20%26%20%3D%20%5Cfrac%7Bg%20%5Calpha%20%5Cbeta%20r_o%5E4%7D%7B%5Cnu%5Ckappa%7D%20%5Cend%7Balign*%7D)

The system is solved using a Toroidal/Poloidal decomposition of the velocity <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\mathbf{u}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{u}" title="\mathbf{u}" /></a>:

<img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{u}=\mathbf{\nabla}\times&space;T&space;\mathbf{r}&space;&plus;&space;\mathbf{\nabla}\times\mathbf{\nabla}\times&space;P&space;\mathbf{r}" title="\mathbf{u}=\mathbf{\nabla}\times T \mathbf{r} + \mathbf{\nabla}\times\mathbf{\nabla}\times P \mathbf{r}" />
![Toroidal/Poloidal decomposition](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cmathbf%7Bu%7D%3D%5Cmathbf%7B%5Cnabla%7D%5Ctimes%20T%20%5Cmathbf%7Br%7D%20&plus;%20%5Cmathbf%7B%5Cnabla%7D%5Ctimes%5Cmathbf%7B%5Cnabla%7D%5Ctimes%20P%20%5Cmathbf%7Br%7D)
