# Double Pendulum


## Introduction
In this project, a double pendulum simulator is implemented using Lagrangian mechanics and a 4th order Runge-Kutta integrator (RK4). The simulation is presented in a graphical interface with real-time animation and plots for angular positions and energy evolution.

**Note: you can find and .exe file on the Release section ready to use!**

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
<img width="341" height="45" alt="image" src="https://github.com/user-attachments/assets/0c2f2a96-17e7-4009-8053-972848691b69" />
</p>
with cartesian coordinates:
<p align="center">
<img width="353" height="75" alt="image" src="https://github.com/user-attachments/assets/5dbd8dbf-1a5e-4db0-84a0-345422ffb72f" />
</p>
The equations of motion are derived using the Euler-Lagrange formalism:
<p align="center">
<img width="341" height="67" alt="image" src="https://github.com/user-attachments/assets/0fcaed2a-ff33-4a68-9106-ff4c6c177979" />
</p>
where b is viscous damping coefficient. The Euler-Lagrange equations lead to a coupled set of nonlinear second-order differential equations:

<p align="center">
<img width="548" height="61" alt="image" src="https://github.com/user-attachments/assets/0fd94a1d-b163-4468-a299-98132676f328" />
</p>

Using a matrix formalism, the angular accelerations can be expressed as:
<p align="center">
<img width="566" height="83" alt="image" src="https://github.com/user-attachments/assets/80d0c715-d816-4c7b-b89b-5285c97b2c5a" />
</p>
where M is the matrix that groups the coefficients of the accelerations, respectively the first two terms of each of the two differential equations written in system form. The matrix C groups the coefficients of the velocities, related to the Coriolis and centripetal effects, while the vector G contains the contributions due to gravity.

Solving for the acceleration, useful for numerical implementation:
<p align="center">
  <img width="159" height="40" alt="image" src="https://github.com/user-attachments/assets/d3cd5360-56a6-4b47-bd56-ec96106afd11" />
</p>

thus:

<p align="center">
<img width="741" height="186" alt="image" src="https://github.com/user-attachments/assets/905460c3-a387-4ddf-9743-02225cd3f969" />
</p>
It is important to note that the explicit form of the accelerations may appear different in some educational texts or books. In fact, the terms are often rewritten by means of algebraic simplifications, collected in different ways or normalised with respect to lengths, masses or numerical constants to obtain more compact or didactic denominators. However, the representation proposed here is formally correct, and above all, directly derivable from the matrix formalism and the previous Euler-Lagrange equations, of
consequence fully verifiable. 

## Numerical integration
For numerical integration, the classical Runge-Kutta 4 method is applied to advance the state vector 𝑦 over discrete time steps Δt (editable in the GUI). The higher the time step, the smoother the animation, but the lower the physical accuracy. Conversely, the lower the time step, the less smooth the animation, but the more physically accurate the simulation. The current time step (0.02) is a good compromise:
<p align="center">
<img width="352" height="225" alt="image" src="https://github.com/user-attachments/assets/5a6d2003-2d3c-4fbc-8fe0-1e47b1f11095" />
</p>
This method ensures a reliable balance between computational efficiency and accuracy, especially in the presence of the system's nonlinear and chaotic behavior. 

## Graphical Interface
The simulator is implemented with Tkinter and features a responsive window that adapts to different screen sizes. The interface is divided into two main regions: a control panel and a visualization panel. The control panel allows adjustment of all physical parameters, initial conditions, and simulation parameters such as time step and friction coefficients. The visualization panel displays the pendulum animation, along with two auxiliary plots: angles versus time and energy versus time.

<img width="872" height="1016" alt="2" src="https://github.com/user-attachments/assets/7ea62775-1ad8-4961-914e-2c0daf9a8cb7" />

The pendulum animation represents each bob as a colored circle connected by rods, with optional trails to visualize the trajectory. The trail lines are semi-transparent and retain a fixed number of past positions, providing insight into the system's chaotic motion. The radius of each bob is proportional to the square root of its mass, providing a direct visual cue for relative mass distribution.

The angles versus time plot wraps angular positions in the [−π,π] range to avoid discontinuities and to facilitate comparison with typical theoretical analyses. The energy plot dynamically scales the vertical axis and displays kinetic, potential, and total energy using distinct colors. These plots are continuously updated alongside the animation, offering a comprehensive view of both the geometrical evolution and the energy dynamics of the system.

## Simulation
The simulation runs in discrete time steps according to the RK4 integrator. The user can interrupt or continue the simulation at any time, reset the system to the initial conditions or switch the visibility of the trajectories. Changes to parameters or initial conditions are made in real time and the simulation is automatically reset to reflect these updates. This interactivity allows experimentation with different masses, lengths, damping and initial angular displacements and provides instant visual feedback on the resulting dynamics.

The combination of real-time animation, angle diagrams and energy monitoring allows a deep understanding of the complex behaviour of the double pendulum, including periodic, quasi-periodic and chaotic states. By slightly changing the initial conditions, users can observe drastically different results, illustrating the sensitive dependence on initial conditions that characterises chaotic systems.

**This simulator represents a bridge between theoretical physics, applied mathematics, and computational implementation, making it suitable for educational, research, and exploratory purposes.**

## Contributions
If you find any bugs, errors or anything else, feel free to open an issue!
