Stick-slip behavior is a phenomenon that occurs when two surfaces slide over each other, and can be described as surfaces alternating between sticking to each other and sliding over each other, with a corresponding change in the force of friction. The static friction coefficient between two surfaces is typically larger than the kinetic friction coefficient, and if an applied force is large enough to overcome the static friction, then the reduction of the friction to the kinetic friction can cause a sudden jump in the velocity, resulting in the jerking motion

In this Julia code we will model a spring mass damper system forced by a sinusoidal forcing term on a rough ground with friction coefficients $ \mu_s $ and $ \mu_k $ being static and kinematic friction coefficients respectively. The equation of motion of the system can then be written as 

$$ m\ddot{x} + c \dot{x} + kx + f_{fric} = f_{ext}  $$

where 

$$
 \ f_{fric} =
  \begin{cases}
      \mu_k N* sgn(\dot{x})     & \quad \text{if } \dot{x} > 0 \\
   min(\mu_s N , f_{ext} - kx)  & \quad \text{if } \dot{x} < 0
  \end{cases}
\
$$


