# ------------------------- #
# DOUBLE PENDULUM SIMULATOR #
# ------------------------- #
import numpy as np
import tkinter as tk    
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.animation import FuncAnimation
from collections import deque

# ---------------------------------------------------------------------------- #
# 1. PHYSICAL-MATHEMATICAL MODEL
# ---------------------------------------------------------------------------- #
class PhysicsModel:
    """
    Contains the physical and numerical model.
    Implements the equations of motion for a double pendulum using Lagrangian mechanics and integrates them using RK4.
    Functions:
        - __init__          -> Physical parameters (m, l, g, b)
        - get_positions     -> Returns the current position of the pendulum
        - compute_alphas    -> Computes angular accelerations from Lagrangian dynamics
        - rk4_step          -> RK4 integrator
        - compute_energies  -> Calculation of positions and energies
    """

    def __init__(self, m1=1.0, l1=1.0, m2=1.0, l2=1.0, g=9.81, b1=0.0, b2=0.0):
        self.m1, self.l1 = m1, l1   # First bob:  mass1 [kg], length1 [m]
        self.m2, self.l2 = m2, l2   # Second bob: mass2 [kg], length2 [m]
        self.g = g                  # g = 9.81 [m/s^2]
        self.b1, self.b2 = b1, b2   # Friction coefficients b1 [N·s/m], b2 [N·s/m]

    def get_positions(self, state):
        """ Calculate the cartesian positions using the state """
        theta1, _, theta2, _ = state

        # Convert angles (theta1, theta2) to cartesian coordinates (x1, y1, x2, y2)
        x1 = self.l1 * np.sin(theta1)       # l * sin(theta1)
        y1 = -self.l1 * np.cos(theta1)      # -l * cos(theta1)
        x2 = x1 + self.l2 * np.sin(theta2)  # x1 + l * sin(theta2)
        y2 = y1 - self.l2 * np.cos(theta2)  # y1 - l * cos(theta2)

        return (x1, y1), (x2, y2)

    def compute_alphas(self, theta1, theta2, omega1, omega2):
        """ Calculate the angular accelerations (domega1/dt, domega2/dt) """
        delta = theta2 - theta1
        cos_delta, sin_delta = np.cos(delta), np.sin(delta)

        den1 = (self.m1 + self.m2) * self.l1 - self.m2 * self.l1 * cos_delta**2
        den2 = (self.l2 / self.l1) * den1
        if abs(den1) < 1e-12: den1 = 1e-12
        if abs(den2) < 1e-12: den2 = 1e-12

        # Numerator for d(omega1)/dt, with addition of the friction term
        domega1_num = (self.m2 * self.l1 * omega1**2 * sin_delta * cos_delta +
                       self.m2 * self.g * np.sin(theta2) * cos_delta +
                       self.m2 * self.l2 * omega2**2 * sin_delta -
                       (self.m1 + self.m2) * self.g * np.sin(theta1) -
                       self.b1 * omega1)  # Viscous damping b1 [N*m*s/rad]

        # Numerator for d(omega2)/dt, with addition of the friction term
        domega2_num = (-self.m2 * self.l2 * omega2**2 * sin_delta * cos_delta +
                       (self.m1 + self.m2) * (self.g * np.sin(theta1) * cos_delta -
                                              self.l1 * omega1**2 * sin_delta -
                                              self.g * np.sin(theta2)) -
                       self.b2 * omega2)  # Viscous damping b2 [N*m*s/rad]

        return domega1_num / den1, domega2_num / den2

    def rk4_step(self, state, dt):
        """ Advance the simulation by a time step 'dt' using Runge-Kutta 4 """
        def f(y):
            theta1, omega1, theta2, omega2 = y
            alpha1, alpha2 = self.compute_alphas(theta1, theta2, omega1, omega2)
            return np.array([omega1, alpha1, omega2, alpha2])

        y = state

        # Compute k1, k2, k3, k4 steps
        k1 = dt * f(y)        
        k2 = dt * f(y + 0.5 * k1) 
        k3 = dt * f(y + 0.5 * k2)  
        k4 = dt * f(y + k3)  

        # Standard RK4 integration y_{n+1} = y_n + (k1 + 2*k2 + 2*k3 + k4)/6 
        return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6  # Next state

    def compute_energies(self, state):
        """ Calculate the kinetic, potential and total energies """
        theta1, omega1, theta2, omega2 = state

        # Kinetic Energy T
        T1 = 0.5 * self.m1 * (self.l1 * omega1)**2  
        T2 = 0.5 * self.m2 * ((self.l1 * omega1)**2 + (self.l2 * omega2)**2 +
                              2 * self.l1 * self.l2 * omega1 * omega2 * np.cos(theta1 - theta2))  

        kinetic = T1 + T2   

        y1 = -self.l1 * np.cos(theta1)      # y1 = -l1 * cos(theta1)
        y2 = y1 - self.l2 * np.cos(theta2)  # y2 = y1 - l2 * cos(theta2)

        # Potential Energy V
        # In most textbooks the formula for V is equal to -(m1 + m2)*g*l1*cos(theta1) - m2*g*l2*cos(theta2), where V=0 is taken at the pivot point.
        # In this simulator, the potential energy V is calculated with respect to the lowest possible point of each bob's swing (where V = 0). This ensures V is always non-negative.
        # Consequently we should add l1 and l2 to the classical potential energy formula
        potential = self.m1 * self.g * (y1 + self.l1) + self.m2 * self.g * (y2 + self.l1 + self.l2)

        return kinetic, potential, kinetic + potential

# ---------------------------------------------------------------------------- #
# 2. GUI AND SIMULATOR
# ---------------------------------------------------------------------------- #
class DoublePendulumSimulator:
    """
    Manages the entire application: the GUI, controls, and animation loop.
    Functions:
        - __init__              -> Initialize simulation, state and GUI elements
        - wrap_angle
        - setup_gui             -> Initialize simulation and main GUI elements
        - setup_controls        -> Create and layout control elements
        - create_param_control  -> Create labeled scale and entry widgets for parameters
        - setup_visualization   -> Sets up the matplotlib figure and subplots
        - setup_plot_elements   -> Initializes all plot lines and axes settings
        - animation_update      -> Update simulation state 
        - apply_changes         -> Apply changes from GUI and resets simulation
        - reset_simulation      -> Reset state, time, and all data arrays
        - toggle_pause          -> Pause / resume the simulation
        - toggle_trails         -> Show or hide the trails
        - run                   -> Start the simulation
    """
    
    def __init__(self):
        # Physics model
        self.physics = PhysicsModel()

        # Initial [theta1, omega1, theta2, omega2]
        self.initial_state = np.array([np.pi / 2, 0.0, np.pi / 2, 0.0]) 
        self.state = self.initial_state.copy()
        self.dt = 0.02                                               # Time step [s]
        self.time = 0.0                                              # Simulation time [s]
        self.running = True                                          # Running flag
        self.show_trails = True                                      # Display trails flag

        # Trail visualization
        self.max_trail_length = 500                                  # Maximum length of the trail
        self.trail1 = deque(maxlen=self.max_trail_length)            # First bob trail
        self.trail2 = deque(maxlen=self.max_trail_length)            # Second bob trail

        # Data storage for plotting
        self.max_data_points = int(20 / self.dt)                     # Max data points for plots
        self.time_data = deque(maxlen=self.max_data_points)          # Time series
        self.theta1_data = deque(maxlen=self.max_data_points)        # Angle theta1
        self.theta2_data = deque(maxlen=self.max_data_points)        # Angle theta2
        self.kinetic_data = deque(maxlen=self.max_data_points)       # Kinetic energy
        self.potential_data = deque(maxlen=self.max_data_points)     # Potential energy
        self.total_data = deque(maxlen=self.max_data_points)         # Total energy

        # Initialize GUI and apply changes
        self.setup_gui()
        self.apply_changes()

    def wrap_angle(self, angle):
        """ Wrap angle to [-π, π] range """
        while angle > np.pi:
            angle -= 2 * np.pi
        while angle < -np.pi:
            angle += 2 * np.pi
        return angle

    def setup_gui(self):
        """ Configure the graphical interface """
        self.root = tk.Tk()
        self.root.title("Double Pendulum Simulation")
        self.root.state("zoomed") 

        # Main frame
        main_frame = ttk.Frame(self.root) 
        main_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Control frame
        control_frame = ttk.LabelFrame(main_frame, text="Controls", padding=10)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 5))

        self.setup_controls(control_frame)  

        # Visualization frame
        viz_frame = ttk.Frame(main_frame)
        viz_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        self.setup_visualization(viz_frame) 

    def setup_controls(self, parent):
        """ Configure the control elements of the interface. """
        row = 0 

        # Simulation controls
        sim_frame = ttk.LabelFrame(parent, text="Simulation", padding=5)
        sim_frame.grid(row=row, column=0, sticky="ew", pady=(0, 10))
        row += 1        

        # Pause / Resume button
        self.pause_btn = ttk.Button(sim_frame, text="Pause", command=self.toggle_pause)
        self.pause_btn.pack(fill=tk.X, pady=2)

        # Reset button
        ttk.Button(sim_frame, text="Reset", command=self.reset_simulation).pack(fill=tk.X, pady=2)

        # Toggle trail visualization
        self.trail_var = tk.BooleanVar(value=self.show_trails)
        ttk.Checkbutton(sim_frame, text="Show Trails", variable=self.trail_var, command=self.toggle_trails).pack(pady=2)

        # Physical parameters
        params_frame = ttk.LabelFrame(parent, text="Physical Parameters", padding=5)
        params_frame.grid(row=row, column=0, sticky="ew", pady=(0, 10))
        row += 1

        # Tkinter variables for physical parameters
        self.m1_var = tk.DoubleVar(value=1.0)
        self.m2_var = tk.DoubleVar(value=1.0)
        self.l1_var = tk.DoubleVar(value=1.0)
        self.l2_var = tk.DoubleVar(value=1.0)
        self.g_var = tk.DoubleVar(value=9.81)
        self.b1_var = tk.DoubleVar(value=0.0)
        self.b2_var = tk.DoubleVar(value=0.0)

        # Create sliders and entry fields for physical parameters
        self.create_param_control(params_frame, "Mass 1 [kg]:", self.m1_var, 0.1, 5.0)
        self.create_param_control(params_frame, "Mass 2 [kg]:", self.m2_var, 0.1, 5.0)
        self.create_param_control(params_frame, "Length 1 [m]:", self.l1_var, 0.1, 2.0)
        self.create_param_control(params_frame, "Length 2 [m]:", self.l2_var, 0.1, 2.0)
        self.create_param_control(params_frame, "Gravity [m/s²]:", self.g_var, 0.1, 20.0)
        self.create_param_control(params_frame, "Friction 1:", self.b1_var, 0.0, 1.0)
        self.create_param_control(params_frame, "Friction 2:", self.b2_var, 0.0, 1.0)

        # Initial conditions frame
        initial_frame = ttk.LabelFrame(parent, text="Initial Conditions", padding=5)
        initial_frame.grid(row=row, column=0, sticky="ew", pady=(0, 10))
        row += 1

        # Tkinter variables for angles and angular velocities
        self.theta1_var = tk.DoubleVar(value=np.degrees(self.initial_state[0]))
        self.theta2_var = tk.DoubleVar(value=np.degrees(self.initial_state[2]))
        self.omega1_var = tk.DoubleVar(value=np.degrees(self.initial_state[1]))
        self.omega2_var = tk.DoubleVar(value=np.degrees(self.initial_state[3]))

        # Create sliders and entry fields for initial conditions
        self.create_param_control(initial_frame, "θ₁ initial [°]:", self.theta1_var, -180, 180)
        self.create_param_control(initial_frame, "θ₂ initial [°]:", self.theta2_var, -180, 180)
        self.create_param_control(initial_frame, "ω₁ initial [°/s]:", self.omega1_var, -360, 360)
        self.create_param_control(initial_frame, "ω₂ initial [°/s]:", self.omega2_var, -360, 360)

        # Simulation parameters frame
        sim_params_frame = ttk.LabelFrame(parent, text="Simulation Parameters", padding=5)
        sim_params_frame.grid(row=row, column=0, sticky="ew", pady=(0, 10))
        row += 1

        # Time step (dt) control
        self.dt_var = tk.DoubleVar(value=self.dt)
        self.create_param_control(sim_params_frame, "Time Step dt [s]:", self.dt_var, 0.001, 0.04)

       # Apply changes button (recreate physics model and reset)
        ttk.Button(parent, text="Apply and Reset", command=self.apply_changes).grid(row=row, column=0, sticky="ew", pady=10)

    def create_param_control(self, parent, label, var, min_val, max_val):
        """ Create a labeled control with sliders and entry fields """

        # Main frame for the control
        frame = ttk.Frame(parent)
        frame.pack(fill=tk.X, pady=2)

        # Label
        ttk.Label(frame, text=label).pack(anchor=tk.W)
        
        # Sub-frame for the slider and entry
        sub_frame = ttk.Frame(frame)
        sub_frame.pack(fill=tk.X)
        
        # Horizontal slider (scale) for adjusting the parameter        
        scale = ttk.Scale(sub_frame, from_=min_val, to=max_val, variable=var, orient=tk.HORIZONTAL)
        scale.pack(fill=tk.X, side=tk.LEFT, expand=True)
        
        # Entry field for precise parameter value
        entry = ttk.Entry(sub_frame, textvariable=var, width=8)
        entry.pack(side=tk.RIGHT, padx=(5, 0))

    def setup_visualization(self, parent):
        """ Set up the graphical visualization with subplot2grid """
        self.fig = plt.figure(figsize=(12, 8))
        
        self.ax_pendulum = plt.subplot2grid((2, 2), (0, 0), rowspan=2)  # Pendulum plot
        self.ax_angles = plt.subplot2grid((2, 2), (0, 1))               # Angles vs time plot
        self.ax_energy = plt.subplot2grid((2, 2), (1, 1))               # Energy vs time plot
        plt.tight_layout(pad=4.0)

        self.canvas = FigureCanvasTkAgg(self.fig, parent)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.setup_plot_elements()

    def setup_plot_elements(self):
        """ Initialize plot elements """

        # Pendulum plot
        self.ax_pendulum.set_title('Double Pendulum')
        self.ax_pendulum.grid(True, alpha=0.3)

        # --- FIX: Separate Rods and Bobs ---
        # 1. Rods (Lines only, black)
        self.rod1, = self.ax_pendulum.plot([], [], '-', linewidth=2, color='black')
        self.rod2, = self.ax_pendulum.plot([], [], '-', linewidth=2, color='black')

        # 2. Bobs (Markers only, colored)
        # zorder=3 ensures bobs are drawn ON TOP of the rods
        self.bob1, = self.ax_pendulum.plot([], [], 'o', color='blue', zorder=3)
        self.bob2, = self.ax_pendulum.plot([], [], 'o', color='red', zorder=3)
        # -----------------------------------

        # Trail lines
        self.trail1_line, = self.ax_pendulum.plot([], [], '-', alpha=0.6, color='blue')
        self.trail2_line, = self.ax_pendulum.plot([], [], '-', alpha=0.6, color='red')

        self.ax_pendulum.set_aspect('equal', adjustable='box')

        # Angles plot
        self.ax_angles.set_title('Angles vs Time (Wrapped [-π, π])')
        self.ax_angles.set_ylabel('Angle (radians)')
        self.ax_angles.grid(True, alpha=0.3)
        self.ax_angles.set_ylim(-np.pi - 0.2, np.pi + 0.2)
        self.theta1_line, = self.ax_angles.plot([], [], label='θ₁', color='blue')
        self.theta2_line, = self.ax_angles.plot([], [], label='θ₂', color='red')
        self.ax_angles.legend()

        # Energy plot
        self.ax_energy.set_title('Energy vs Time')
        self.ax_energy.set_xlabel('Time (s)')
        self.ax_energy.set_ylabel('Energy (J)')
        self.ax_energy.grid(True, alpha=0.3)
        self.kinetic_line, = self.ax_energy.plot([], [], label='Kinetic', color='green')
        self.potential_line, = self.ax_energy.plot([], [], label='Potential', color='orange')
        self.total_line, = self.ax_energy.plot([], [], label='Total', color='purple', lw=2)
        self.ax_energy.legend()

    def animation_update(self, frame):
        """ Update function called by FuncAnimation for each frame """
        
        # Step simulation if running
        if self.running:
            self.state = self.physics.rk4_step(self.state, self.dt) # RK4 integration step
            self.time += self.dt                                    # Increment simulation time by dt

        (x1, y1), (x2, y2) = self.physics.get_positions(self.state)
        
        # Update Rods
        self.rod1.set_data([0, x1], [0, y1])
        self.rod2.set_data([x1, x2], [y1, y2])

        # Update Bobs (Size depends on their mass)
        self.bob1.set_data([x1], [y1])
        self.bob1.set_markersize(8 * np.sqrt(self.physics.m1)) 

        self.bob2.set_data([x2], [y2])
        self.bob2.set_markersize(8 * np.sqrt(self.physics.m2)) 
        # --------------------------------------------

        # Update trails if running and enabled
        if self.running and self.show_trails:
            # Append positions to trails
            self.trail1.append((x1, y1))
            self.trail2.append((x2, y2))

            # Update trail line if enough length
            if len(self.trail1) > 1: 
                self.trail1_line.set_data(*zip(*self.trail1))
            if len(self.trail2) > 1: 
                self.trail2_line.set_data(*zip(*self.trail2))

        # Update angles and energy data
        if self.running:
            # Wrap angles to [-pi, pi] range for display
            wrapped_theta1 = self.wrap_angle(self.state[0])
            wrapped_theta2 = self.wrap_angle(self.state[2])
            
            # Append time series and wrapped angular positions data
            self.time_data.append(self.time)                        
            self.theta1_data.append(wrapped_theta1)                   
            self.theta2_data.append(wrapped_theta2)                   
            
            ke, pe, te = self.physics.compute_energies(self.state)  
            
            # Append energies
            self.kinetic_data.append(ke)
            self.potential_data.append(pe)
            self.total_data.append(te)

            # Update plots lines
            self.theta1_line.set_data(self.time_data, self.theta1_data)
            self.theta2_line.set_data(self.time_data, self.theta2_data)
            self.kinetic_line.set_data(self.time_data, self.kinetic_data)
            self.potential_line.set_data(self.time_data, self.potential_data)
            self.total_line.set_data(self.time_data, self.total_data)

            # Scale axes (but not angles plot, which has fixed limits)
            for ax in [self.ax_energy]:
                ax.relim(); ax.autoscale_view()
                
            # Scale only x-axis for angles plot to show time progression
            if len(self.time_data) > 1:
                self.ax_angles.set_xlim(min(self.time_data), max(self.time_data))
        
        # Return all artists to re-draw
        return (self.rod1, self.rod2, self.bob1, self.bob2, 
                self.trail1_line, self.trail2_line,
                self.theta1_line, self.theta2_line, self.kinetic_line,
                self.potential_line, self.total_line)

    def apply_changes(self):
        """ Apply changes from the control elements and reset the simulation """
        
        # Update dt from the slider
        self.dt = self.dt_var.get()

        # Recreate the physical model with the updated parameters
        self.physics = PhysicsModel(
            m1=self.m1_var.get(), l1=self.l1_var.get(), # mass and length of first pendulum
            m2=self.m2_var.get(), l2=self.l2_var.get(), # mass and length of second pendulum   
            g=self.g_var.get(),                         # g
            b1=self.b1_var.get(), b2=self.b2_var.get()  # damping coefficients
        )

        # Set initial state
        self.initial_state = np.array([
            np.radians(self.theta1_var.get()), np.radians(self.omega1_var.get()),
            np.radians(self.theta2_var.get()), np.radians(self.omega2_var.get())])

        # Update max data points based on new dt
        self.max_data_points = int(20 / self.dt)

        # Set axis limits
        limit = (self.physics.l1 + self.physics.l2) * 1.1
        self.ax_pendulum.set_xlim(-limit, limit)
        self.ax_pendulum.set_ylim(-limit, limit)
        self.reset_simulation()

    def reset_simulation(self):
        """ Reset the simulation to initial conditions """
        self.state = self.initial_state.copy()
        self.time = 0.0

        # Recreate deques with updated max length
        self.time_data = deque(maxlen=self.max_data_points)
        self.theta1_data = deque(maxlen=self.max_data_points)
        self.theta2_data = deque(maxlen=self.max_data_points)
        self.kinetic_data = deque(maxlen=self.max_data_points)
        self.potential_data = deque(maxlen=self.max_data_points)
        self.total_data = deque(maxlen=self.max_data_points)

        # Clear trails
        self.trail1.clear()
        self.trail2.clear()
        
        # Clear plot lines
        self.trail1_line.set_data([], [])
        self.trail2_line.set_data([], [])

        if not self.running: self.toggle_pause()

    def toggle_pause(self):
        """ Toggle the running state """
        self.running = not self.running
        self.pause_btn.config(text="Resume" if not self.running else "Pause")
        
    def toggle_trails(self):
        """ Toggle trails """
        self.show_trails = self.trail_var.get()

        # If trails are not shown, clear the trail data
        if not self.show_trails:
            self.trail1.clear(); self.trail2.clear()
            self.trail1_line.set_data([], []); self.trail2_line.set_data([], [])

    def run(self):
        """ Run the simulation and animation """
        self.animation = FuncAnimation(
            self.fig,
            self.animation_update,
            interval=int(self.dt*1000),
            blit=True,
            cache_frame_data=False
        )
        self.root.mainloop()

if __name__ == "__main__":
    simulator = DoublePendulumSimulator()
    simulator.run()
