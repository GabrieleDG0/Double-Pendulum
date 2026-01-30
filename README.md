# Double Pendulum - Lagrangian Mechanics RK4 Integration

**Note: you can find and .exe file on the Release section ready to use!**
**Note: While the project is fully functional, the current code requires refactoring**

## Introduction
In this project, a double pendulum simulator is implemented using Lagrangian mechanics and a 4th order Runge-Kutta integrator (RK4). The simulation is presented in a graphical interface with real-time animation and plots for angular positions and energy evolution.

<img width="1919" height="1016" alt="dp-simulator" src="https://github.com/user-attachments/assets/2157765e-9cce-4864-bd38-75ee9c911e18" />

## Physical and Mathematical model
The double pendulum consists of two rigid rods, each with a concentrated mass at its end, connected serially with the first rod pivoted to a fixed point. The system is governed by Lagrangian mechanics, where the Lagrangian L is the difference between the kinetic energy T and the potential energy V:
<p align="center">
  <img width="135" height="50" alt="image" src="https://github.com/user-attachments/assets/47f0d3aa-4b22-49db-a7f2-6b8d03175150" />
</p>
The kinetic energy accounts for the translational motion of both bobs, including the interaction between the two masses:
<p align="center">
<img width="592" height="67" alt="image" src="https://github.com/user-attachments/assets/2a906099-b9c0-41a9-baa0-ca460d7f2b1c" />
</p>
The potential energy is determined by the height of each mass in the gravitational field:
<p align="center">
<img width="392" height="51" alt="image" src="https://github.com/user-attachments/assets/8fb68a24-f049-47ef-9744-68c073990738" />
</p>
with cartesian coordinates:
<p align="center">
<img width="353" height="75" alt="image" src="https://github.com/user-attachments/assets/5dbd8dbf-1a5e-4db0-84a0-345422ffb72f" />
</p>
The Lagrangian of the system L = T - V:
<p align="center">
<img width="560" height="97" alt="image" src="https://github.com/user-attachments/assets/9628c716-dbf3-4975-8b6e-97c961d5c04b" />
</p>
The equations of motion are derived using the Euler-Lagrange formalism:
<p align="center">
<img width="341" height="67" alt="image" src="https://github.com/user-attachments/assets/0fcaed2a-ff33-4a68-9106-ff4c6c177979" />
</p>
where b is viscous damping coefficient. The Euler-Lagrange equations lead to a coupled set of nonlinear second-order differential equations:

<p align="center">
<img width="698" height="87" alt="image" src="https://github.com/user-attachments/assets/df7f2588-b1a5-41ae-a375-16ed816d904b" />
</p>

Using a matrix formalism, the angular accelerations can be expressed as:
<p align="center">
<img width="789" height="119" alt="image" src="https://github.com/user-attachments/assets/fe71f428-1451-482b-a14e-5bf3292f1c96" />
</p>
where M is the matrix that groups the coefficients of the accelerations, respectively the first two terms of each of the two differential equations written in system form. The matrix C groups the coefficients of the velocities, related to the Coriolis and centripetal effects, while the vector G contains the contributions due to gravity.

Solving for the acceleration, useful for numerical implementation:
<p align="center">
  <img width="206" height="55" alt="image" src="https://github.com/user-attachments/assets/c6e3adca-fbad-4765-82cf-0b65959bfa88" />
</p>

thus:

<p align="center">
<img width="933" height="238" alt="image" src="https://github.com/user-attachments/assets/0e867252-2b19-4e7f-97d9-c1d757a102f0" />
</p>
It is important to note that the explicit form of the accelerations may appear different in some educational texts or books. In fact, the terms are often rewritten by means of algebraic simplifications, collected in different ways or normalised with respect to lengths, masses or numerical constants to obtain more compact or didactic denominators. However, the representation proposed here is formally correct, and above all, directly derivable from the matrix formalism and the previous Euler-Lagrange equations, of
consequence fully verifiable. 

## Numerical integration
The 4th-order Runge-Kutta is a numerical further, in particular a four-stage explicit integration method that combines several evaluations of the f-function in the step to achieve higher
precision. Four increments (or slopes) are calculated at each interval h: 
<p align="center">
<img width="288" height="189" alt="image" src="https://github.com/user-attachments/assets/7b01dc0a-98c9-46bf-bbaa-b054022d7c59" />
</p>
The next step approximation is obtained as a weighted combination of these values:
<p align="center">
<img width="356" height="66" alt="image" src="https://github.com/user-attachments/assets/64d7f441-07c4-47f8-aa23-a067525c727c" />
</p>
The higher the time step, the smoother the animation, but the lower the physical accuracy. Conversely, the lower the time step, the less smooth the animation, but the more physically accurate the simulation. The current time step (0.02) is a good compromise.
In practical terms, RK4 achieves high accuracy with a much smaller number of steps than Euler method. The disadvantage is that it requires four times as many evaluations of the function
f per step (4 steps instead of 1), but on the whole this additional cost is compensated for by the possibility of using considerably larger h. The RK4 method is in fact considered an excellent compromise between accuracy, computational cost and stability, even if it is more expensive per step. 
Graphical RK4 on approximating the blue function:
<p align="center">
<img width="453" height="460" alt="image" src="https://github.com/user-attachments/assets/5d5ae645-3833-44f9-8029-556793a65177" />
</p>

## Graphical Interface
The simulator is implemented with Tkinter and features a responsive window that adapts to different screen sizes. The interface is divided into two main regions: a control panel and a visualization panel. The control panel allows adjustment of all physical parameters, initial conditions, and simulation parameters such as time step, the friction coefficients, masses, lenghts and angles. The visualization panel displays the pendulum animation, along with two auxiliary plots: angles versus time and energy versus time.

<img width="872" height="1016" alt="2" src="https://github.com/user-attachments/assets/7ea62775-1ad8-4961-914e-2c0daf9a8cb7" />

The pendulum animation represents each bob as a colored circle connected by rods, with optional trails to visualize the trajectory. The trail lines are semi-transparent and retain a fixed number of past positions, providing insight into the system's chaotic motion. The radius of each bob is proportional to the square root of its mass, providing a direct visual cue for relative mass distribution.

The angles versus time plot wraps angular positions in the [−π,π] range to avoid discontinuities and to facilitate comparison with typical theoretical analyses. The energy plot dynamically scales the vertical axis and displays kinetic, potential, and total energy using distinct colors. These plots are continuously updated alongside the animation, offering a comprehensive view of both the geometrical evolution and the energy dynamics of the system.

## Simulation
The simulation runs in discrete time steps according to the RK4 integrator. The user can interrupt or continue the simulation at any time, reset the system to the initial conditions or switch the visibility of the trajectories. Changes to parameters or initial conditions are made in real time and the simulation is automatically reset to reflect these updates. This interactivity allows experimentation with different masses, lengths, damping and initial angular displacements and provides instant visual feedback on the resulting dynamics.

The combination of real-time animation, angle diagrams and energy monitoring allows a deep understanding of the complex behaviour of the double pendulum, including periodic, quasi-periodic and chaotic states. By slightly changing the initial conditions, users can observe drastically different results, illustrating the sensitive dependence on initial conditions that characterises chaotic systems.

**This simulator represents a bridge between theoretical physics, applied mathematics, and computational implementation, making it suitable for educational, research, and exploratory purposes.**

### Contributions
If you find any bugs, errors or anything else, feel free to open an issue!
