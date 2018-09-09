# BBDYN
Repository includes source code for a mechanistic agent-based bark beetle model BBDYN included in WINDROT framework

The code is written in Fortran and includes subroutines and functions necessary for BBDYN model. 
The main subroutine is called in this article by the WINDROT framework. 

N.B. The source code itself is not a functioning program.

The model itself is described in detail in an scientific article:
Honkaniemi, J., Ojansuu, R., Kasanen, R., & Heliövaara, K. 2018. Interaction of disturbance agents on Norway spruce: a mechanistic model of bark beetle dynamics integrated in simulation framework WINDROT. in press Ecological Modelling.

BBDYN_module.for includes

subroutine bbdyn : the main subroutine that is called by the WINDROT framework connecting the forest dynamics and disturbance models together
subroutine bbredistri : redistributes the bark beetles from their initial landing coordinates to nearby trees
subroutine bb_reproduction : Colonization of trees and bark beetle reproduction in succesfully colonized trees 
subroutine indexx1 : Indexes an array in a specified order    
    
function defth : Defence threshold based on tree vigour index 
function Abark : Calculates suitable bark area for bark beetles for each tree
function diarel : Calculates diameter of tree at height l based on stem taper curve by Laasasenaho (1982)
function antag_curve : Probability curve for the antagonist effect
function dheart : Calculates the heartwood diameter at diameter d
