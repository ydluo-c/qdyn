"""
                    .: PyQDYN: Python wrapper for QDYN :.

Author: Martijn van den Ende

This wrapper can be imported and called from a different Python script file to
do the heavy lifting regarding the generation and population of the mesh,
writing of the qdyn input file, and reading of the qdyn output files. See the
examples directory for usage.
"""

import os
import re
import subprocess

import numpy as np
from scipy.optimize import root
import pandas as pd
from pandas import read_csv


# import antigravity	# xkcd.com/353/

__version__ = "3.0.0"


class qdyn:

    mesh_rendered = False
    qdyn_input_written = False
    mesh_dict = {}

    # Directory path to qdyn executable
    # Default: same directory as this python script
    qdyn_path = os.path.abspath(
        os.path.join(os.path.realpath(__file__), os.pardir)
    )

    # Working directory can be kept empty, except for special cases
    work_dir = ""
    # Flag for using the bash environment in Windows 10
    W10_bash = False

    # Initialise qdyn wrapper
    def __init__(self):
        set_dict = {}

        set_dict["FRICTION_MODEL"] = "RSF"	# Model for fault rheology: "RSF" = rate-and-state friction, "CNS" = Chen-Niemeijer-Spiers microphysical model

        # Fault properties and features
        set_dict["MESHDIM"] = 0				# Dimensionality of fault: 0 = spring-slider, 1 = 1D fault in 2D medium, 2 = 2D fault in 3D medium

        # Along strike boundary conditions:
        # 0 = periodic
        # 1 = fixed-length asperity surrounded by steady creep
        # 2 = periodic and symmetric around the first element
        # 3 = fixed-length and symmetric around the first element
        set_dict["FINITE"] = 1

        # ODE solver mode
        # 1 = Bulirsch-Stoer
        # 2 = Runge-Kutta-Fehlberg
        set_dict["SOLVER"] = 1

        # Fault loading geometry
        #  1 = strike-slip (right-lateral)
        # -1 = strike-slip (left-lateral)
        #  2 = thrust
        # -2 = normal
        set_dict["FAULT_TYPE"] = 1

        set_dict["V_PL"] = 1e-9             # Tectonic loading rate [m/s]
        set_dict["L"] = 1					# Length of the fault if MESHDIM = 1. Stiffness is MU/L if MESHDIM = 0
        set_dict["W"] = 50e3	       		# Distance between load point and fault if MESHDIM = 1 and FINITE = 0
        set_dict["SIGMA"] = 50e6            # Normal stress [Pa]
        set_dict["VS"] = 3000				# Shear wave speed [m/s]. If VS = 0 radiation damping is off
        set_dict["MU"] = 30e9				# Shear modulus [Pa]
        set_dict["LAM"] = 30e9				# Elastic modulus for 3D simulations [Pa]
        set_dict["TPER"] = 31536000			# Period of additional time-dependent oscillatory shear stress loading [s]
        set_dict["APER"] = 0				# Amplitude of additional time-dependent oscillatory shear stress loading [Pa]
        set_dict["DIP_W"] = 90				# Fault dip in 3D
        set_dict["Z_CORNER"] = 0			# 3D (Base of the fault)
        set_dict["FAULT_LABEL"] = 1         # Label of the fault (int) 
        set_dict["N_FAULTS"] = -1           # Number of faults (int); this should be set automatically

        # Optional simulation features
        set_dict["FEAT_STRESS_COUPL"] = 0	# Normal stress coupling
        set_dict["FEAT_TP"] = 0 	# Thermal pressurisation
        set_dict["FEAT_FLUID_DIFF"] = 0		# Fluid diffusion
        set_dict["FEAT_VAR_K"] = 0		    # Variable permeability
        set_dict["FEAT_LOCALISATION"] = 0	# Gouge zone localisation of strain (CNS only)
        set_dict["RESTART_STEP"] = 0        # Restart step of the simulation [-]
        set_dict["RESTART_TIME"] = 0        # Restart time of the simulation [s]
        set_dict["RESTART_SLIP"] = 0        # Restart slip of the simulation [m] 

        # Rate-and-state friction parameters
        set_dict["SET_DICT_RSF"] = {
            "RNS_LAW": 0,						# Type of rate-and-state friction law (0: original, 1: with cut-off velocities)
            "THETA_LAW": 1,						# State evolution law (0: aging in no-healing approximation, 1: ageing law, 2: slip law)
            "A": 0.0010,						# Direct effect parameter
            "B": 0.0011,						# Evolution effect parameter
            "DC": 1e-5,							# Characteristic slip distance
            "V1": 0.01,							# Cut-off velocity of direct effect (when RNS_LAW = 1) [m/s]
            "V2": 1e-7,							# Cut-off velocity of evolution effect (when RNS_LAW = 1) [m/s], Controls transition from weakening to strengthening. V2 should be <= V1
            "MU_SS": 0.6,						# Reference steady-state friction
            "V_SS": 1e-6,						# Reference steady-state slip velocity [m/s]
            "CO": 0,							# Cohesion [Pa]
            "TH_0": 1.0,						# Initial state [s]
            "TAU": -1.0,                        # Initial stress [Pa] (default -1 indicating it will be computed later)
            "INV_VISC": 0.0,                    # Inverse viscosity [m/Pa.s]
        }

        # CNS microphysical model parameters
        set_dict["SET_DICT_CNS"] = {
            "THICKNESS": 1e-3,					# Total thickness of gouge zone
            # Granular flow parameters
            "A_TILDE": 0.006,					# Coefficient of logarithmic rate-dependence of grain-boundary friction
            "MU_TILDE_STAR": 0.4,				# Reference grain-boundary friction at Y_GR_STAR
            "Y_GR_STAR": 30e-10,				# Reference strain rate corresponding with MU_TILDE_STAR
            "H": 0.577,							# Dilatancy geometric factor
            # Porosity parameters
            "PHI_C": 0.32,						# Critical state porosity
            "PHI0": 0.03,						# Lower cut-off velocity (e.g. percolation threshold)
            "PHI_INI": 0.25,					# Simulation initial porosity
            # Creep mechanism parameters (grouped in list)
            "N_CREEP": -1,                      # Number of creep mechanisms (must be set by render_mesh() )
            "A": [0.0],				            # (T-dependent) kinetic parameters for creep mechanisms
            "N": [1.0],                         # Stress exponents for creep mechanisms
            "M": [2.0],				            # Porosity exponents for creep mechanisms
            # State of stress
            "TAU": 2e7,							# Initial shear stress [Pa]
        }

        # Localisation settings
        set_dict["SET_DICT_LOCALISATION"] = {
            "LOCALISATION": 1.0,				# Degree of localisation: 1.0 = entire gouge layer is active, 0.0 = entire gouge zone is passive
            "PHI_INI_BULK": 0.2,				# Initial porosity passive gouge zone
            # Creep mechanism parameters for passive/bulk zone (grouped in list)
            "A": [0.0],				            # (T-dependent) kinetic parameters for creep mechanisms
            "N": [1.0],                         # Stress exponents for creep mechanisms
            "M": [2.0],				            # Porosity exponents for creep mechanisms
        }

        # Thermal pressurisation model
        set_dict["SET_DICT_TP"] = {
            "RHOC": 2.16e6,						# Density times specific heat capacity (rho*cp) [J/k/m3]
            "HALFW": 0.5e-3,					# Half-width of the fault zone [m]
            "BETA": 2e-9,						# Bulk compressibility (fluid + pore) [1/Pa]
            "ETA": 2e-4,						# Fluid dynamic viscosity [Pa s]
            "LAM": 2.78e-4,						# Nett thermal expansion coefficient (fluid - solid) [1/K]
            "K_T": 2.0,							# Thermal conductivity [J/s/K/m]
            "K_P": 1e-16,						# Intrinsic hydraulic permeability [m2]
            "P_A": 0.0,							# Ambient fluid pressure [Pa]
            "T_A": 293.0,						# Ambient temperature [K]
            "DILAT_FACTOR": 0.0,                # Factor > 0 to control amount of dilatancy hardening
        }
        #----------------------------------------------------------------------------------
        # Fluid diffusion model 
        set_dict["SET_DICT_FLUID_DIFF"] = {
            "RHOF": 2.16e6,						# Density Fluid
            "BETA": 2e-9,						# Bulk compressibility (fluid + pore) [1/Pa]
            "ETA": 2e-4,						# Fluid dynamic viscosity [Pa s]
            "PHI": 2e-4,                       # Porosity 
            "PERMEABILITY": 1e-14,             # Permeability
            "P_A": 0.0,							# Ambient fluid pressure [Pa]
            "nb_source": 1,                   # Number of source (only one for now)
            "rate_injection":  [1.25e-6],       # Injection rate (m/s)
            "index_injection": [1000],          # Indec of the fault where injection is performed (for FORTRAN NOT PYTHON)
            "t_injection_beg": [0.0],           # Beginning time of injection
            "t_injection_end": [100.0*86400.0]  # End time of injection
        }
        
        set_dict["SET_DICT_VAR_K"] = {
            "Lk": 1.0,
            "Tk": 1000.0,
            "kmin": 1e-13,      
            "kmax": 1e-15,
            "Snk": 1e6,
        }
        #----------------------------------------------------------------------------------

              
                   
                    
                   
        

        # Benjamin Idini's damage model
        set_dict["D"] = 0
        set_dict["HD"] = 0

        # Discretisation and accuracy parameters
        set_dict["N"] = -1					# Number of mesh elements if MESHDIM = 1
        set_dict["NX"] = 1					# Number of fault elements along-strike in 3D
        set_dict["NW"] = -1					# Number of fault elements along-dip in 3D
        set_dict["DW"] = 1					# Length of segment along-strike in 3D
        set_dict["TMAX"] = 1000				# Total simulation time [s]
        set_dict["NSTOP"] = 0				# Stopping criterion (0: stop at t=TMAX, 1: stop end localisation phase, 2: stop first slip rate peak, 3: stop v > TMAX)
        set_dict["DTTRY"] = 1e-1			# First trial step [s]
        set_dict["DTMAX"] = 0				# Maximum time step [s] (0 = unrestricted)
        set_dict["ACC"] = 1e-7				# Solver accuracy

        # Output control parameters
        set_dict["V_TH"] = 1e-2				# Threshold velocity for seismic event
        set_dict["NTOUT_LOG"] = 10          # Temporal interval (number of time steps) for time series output
        set_dict["NTOUT_OT"] = 10			# Temporal interval (number of time steps) for time series output
        set_dict["NTOUT_OX"] = 1000			# Temporal interval (number of time steps) for snapshot output
        set_dict["NXOUT_OX"] = 1			# Spatial interval (number of elements in x-direction) for snapshot output
        set_dict["NWOUT_OX"] = 1            # Spatial interval (number of elements in y-direction) for snapshot output
        set_dict["OX_SEQ"] = 0				# Type of snapshot outputs (0: all snapshots in single file, 1: one file per snapshot)
        set_dict["OX_DYN"] = 0				# Output specific snapshots of dynamic events defined by thresholds on peak slip velocity DYN_TH_ON and DYN_TH_OFF
        set_dict["NXOUT_DYN"] = 1			# Spatial interval (number of elements in x-direction) for dynamic snapshot outputs
        set_dict["NWOUT_DYN"] = 1           # Spatial interval (number of elements in y-direction) for dynamic snapshot outputs
        set_dict["DYN_TH_ON"] = 1e-3		# peak slip rate threshold defining beginning of dynamic event
        set_dict["DYN_TH_OFF"] = 1e-4		# peak slip rate threshold defining end of dynamic event
        set_dict["IC"] = 0					# Index of selected elements for time series output (starting at 0)
        set_dict["IOT"] = 0					# Indices of elements for additional time series outputs
        set_dict["IASP"] = 0				# Flags for elements (identification only)
        set_dict["SUFFIX"] = ""

        # SPECFEM3D parameters
        set_dict["DYN_FLAG"] = 0			# Integration with dynamic code
        set_dict["DYN_SKIP"] = 0			# Number of dynamic events to skip (warm-up cycles)
        set_dict["DYN_M"] = 1e18			# Target seismic moment of dynamic event

        set_dict["NPROC"] = 1				# Number of processors, default 1 = serial (no MPI)
        set_dict["MPI_PATH"] = "mpiexec"     # Path to MPI executable. Default is set to just mpiexec, but this can be an absolute path
        set_dict["VERBOSE"] = 0
        set_dict["DEBUG"] = 0

        self.set_dict = set_dict

        pass

    # Receive and store settings dict given by script file
    def settings(self, set_dict):
        for key in set_dict:
            val = set_dict[key]
            if type(val) == dict:
                self.set_dict[key].update(val)
            else:
                self.set_dict[key] = val
        return True

    # Render the mesh, and populate it uniformly with the input parameters
    # After render_mesh() was called, the self.mesh_dict can be modified
    # in the script file to introduce heterogeneity in the parameters
    # See examples/PyQDYN/Tse+Rice/tse_rice.py for an example
    def render_mesh(self):

        settings = self.set_dict
        dim = settings["MESHDIM"]
        N = settings["N"]

        # Determine number of creep mechanisms
        A = np.array(settings["SET_DICT_CNS"]["A"])
        N_creep = len(A)
        settings["SET_DICT_CNS"]["N_CREEP"] = N_creep

        if dim == 0:
            N = 1
            settings["N"] = N
            settings["NW"] = N

        if dim == 1:
            N = settings["N"]
            settings["NW"] = N

        if dim == 2:
            N = settings["NX"] * settings["NW"]
            settings["N"] = N

        assert N >= 1, "Number of mesh elements needs to be set before rendering mesh, unless MESHDIM = 0 (spring-block)"

        # NOTE: even though TAU is a common parameter independent 
        # of the friction law, it is included in the friction-law 
        # specific dictionaries for compatibility reasons
        mesh_params_general = ("IOT", "IASP", "V_PL", "SIGMA", "TAU")
        mesh_params_RSF = ("TH_0", "A", "B", "DC", "V1", "V2", "MU_SS", "V_SS", "CO", "INV_VISC")
        mesh_params_CNS = ("PHI_INI", "A_TILDE", "MU_TILDE_STAR", "Y_GR_STAR", "H", "PHI_C", "PHI0", "THICKNESS")
        mesh_params_creep = ("A", "N", "M")
        mesh_params_localisation = ("LOCALISATION", "PHI_INI_BULK")
        mesh_params_TP = ("RHOC", "BETA", "ETA", "HALFW", "K_T", "K_P", "LAM", "P_A", "T_A", "DILAT_FACTOR")
        mesh_params_fluid_diff = ( "RHOF","BETA","ETA","PHI","PERMEABILITY", "P_A")
        mesh_params_var_k = ("Tk","Lk","Snk","kmin","kmax")


        mesh_dict = {}

        # Check if the shear stress was set directly to
        # the main settings dictionary. If not: check
        # the friction law dictionaries
        if settings.get("TAU") is None:
            # If RSF is used
            if settings["FRICTION_MODEL"] == "RSF":
                # Add initial velocity as a parameter
                mesh_params_RSF += ("V_0",)
                # Check if the RSF dictionary has TAU specified
                if settings["SET_DICT_RSF"].get("TAU") is None:
                    # If not, check if V_0 is specified (backward compatibility)
                    assert settings["SET_DICT_RSF"].get("V_0") is not None, "Either V_0 or TAU needs to be specified"
                    # Set TAU to -1 to indicate that it needs to be computed later
                    settings["TAU"] = -1
                else:
                    settings["TAU"] = settings["SET_DICT_RSF"]["TAU"]
            
            # If CNS is used
            if settings["FRICTION_MODEL"] == "CNS":
                assert settings["SET_DICT_CNS"].get("TAU"), "TAU has not been specified"
                # Copy TAU from the CNS dictionary
                settings["TAU"] = settings["SET_DICT_CNS"]["TAU"]
            

        # Populate mesh with general settings
        for param in mesh_params_general:
            mesh_dict[param] = np.ones(N)*settings[param]

        # Add time series (IC) output point
        mesh_dict["IOT"][settings["IC"]] = 1

        # Check for friction model
        if settings["FRICTION_MODEL"] == "RSF":
            for param in mesh_params_RSF:
                mesh_dict[param] = np.ones(N)*settings["SET_DICT_RSF"][param]
        elif settings["FRICTION_MODEL"] == "CNS":
            for param in mesh_params_CNS:
                mesh_dict[param] = np.ones(N)*settings["SET_DICT_CNS"][param]
            for i in range(N_creep):
                for param in mesh_params_creep:
                    param_no = "%s%i" % (param, i)
                    mesh_dict[param_no] = np.ones(N)*settings["SET_DICT_CNS"][param][i]
        else:
            raise ValueError("FRICTION_MODEL '%s' not recognised. Supported models: RSF, CNS" % (settings["FRICTION_MODEL"]))

        # Populate mesh with localisation parameters
        if settings["FEAT_LOCALISATION"] == 1:
            for param in mesh_params_localisation:
                mesh_dict[param] = np.ones(N)*settings["SET_DICT_LOCALISATION"][param]
            for i in range(N_creep):
                for param in mesh_params_creep:
                    param_bulk = "%s%i_bulk" % (param, i)
                    mesh_dict[param_bulk] = np.ones(N)*settings["SET_DICT_LOCALISATION"][param][i]

        # Populate mesh with thermal pressurisation parameters
        mesh_dict["P_A"] = np.zeros(N)
        if settings["FEAT_TP"] == 1:
            for param in mesh_params_TP:
                mesh_dict[param] = np.ones(N)*settings["SET_DICT_TP"][param]
                
                
        if settings["FEAT_FLUID_DIFF"] == 1:
            for param in mesh_params_fluid_diff:
                mesh_dict[param] = np.ones(N)*settings["SET_DICT_FLUID_DIFF"][param]
                
        if settings["FEAT_VAR_K"] == 1:
            for param in mesh_params_var_k:
                mesh_dict[param] = np.ones(N)*settings["SET_DICT_VAR_K"][param]
                
                

        # Mesh XYZ coordinates (and dip angle)
        mesh_dict["X"] = np.zeros(N)
        mesh_dict["Y"] = np.zeros(N)
        mesh_dict["Z"] = np.zeros(N)
        mesh_dict["DIP_W"] = np.zeros(N)

        dx = 1.0 * settings["L"] / settings["NX"]
        halfL = settings["L"] / 2.0

        if dim == 1:
            mesh_dict["X"] = np.linspace(-halfL + 0.5*dx, halfL - 0.5*dx, N)
            mesh_dict["Y"] = np.ones(N)*settings["W"]

        if dim == 2:
            # Generate a default mesh with constant dip angle
            # This can be modified with compute_mesh_coords
            dw = settings["W"] / settings["NW"]
            theta = settings["DIP_W"] * np.pi / 180.0
            mesh_dict["DIP_W"] = np.ones(N) * settings["DIP_W"]
            mesh_dict["DW"] = np.ones(N) * dw
            x = (np.arange(settings["NX"]) + 0.5) * dx
            y = (np.arange(settings["NW"]) + 0.5) * dw * np.cos(theta)
            X, Y = np.meshgrid(x, y, indexing="xy")
            mesh_dict["X"] = X.ravel()
            mesh_dict["Y"] = Y.ravel()
            mesh_dict["Z"] = settings["Z_CORNER"] + mesh_dict["Y"]*np.tan(theta)
        
        # Populate mesh with fault labels
        mesh_dict["FAULT_LABEL"] = np.ones(N) * settings["FAULT_LABEL"]

        # Calculate number of faults
        mesh_dict["N_FAULTS"] = len(np.unique(settings["FAULT_LABEL"]))

        # Populate mesh with values of restart slip
        mesh_dict["RESTART_SLIP"] = np.ones(N) * settings["RESTART_SLIP"]

        self.mesh_dict.update(mesh_dict)
        self.mesh_rendered = True

        return True

    def compute_mesh_coords(self, mesh_dict, dip, dw=None):
        """
        This function computes the 3D cartesian coordinates of the discretised
        fault nodes (centres) and overwrites the mesh dictionary. Since the
        dip can be variable (only along-dip), the mesh coordinates are set iteratively
        according to the following procedure:

        1. The mesh is constructed row-by-row with contiguous along-strike entries (FFT-axis)
        2. The base of the fault is located as a depth Z_CORNER
        3. Starting at the base of the fault, the coordinates of each subsequent row are
           computed with a fixed dip angle and along-dip spacing. The along-dip spacing does
           not need to be uniform or continuous

        IMPORTANT: the first element of the dip angle and spacing vectors (`dip` and `dw`)
        are taken to be at the base of the fault, located at Z_CORNER depth. Moreover, the 
        use of the FFT requires that the along-strike direction is uniform and co-linear, 
        which implies that the strike be constant (zero) and that the dip be constant 
        along-strike. Hence, an initial check is performed to verify that the length of the 
        input vectors equals zero (scalar value) or equals Nw (number of elements along-dip).

        Adopted convention for the coordinate frame:
          x = along-strike direction
          y = perpendicular to x in the horizontal plane
          z = vertical (negative from the surface down)
          (w = along-dip direction)

        Adopted convention for the indexation:
          matrix[i, j] = vector[i*Nx + j]  (0 <= i < Nw, 0 <= j < Nx)
          vector[n] = matrix[n // Nx, n % Nx]  (0 <= n < Nw * Nx)
          matrix.shape = (Nw, Nx)
          vector.shape = (Nw * Nx,)
          Nw: elements along-dip
          Nx: elements along-strike

        Starting point of the mesh (first element): 
        (x, y, z) = ( dx/2, dw[0]*cos(dip[0])/2, Z_CORNER + dw[0]*sin(dip[0])/2 )
        """

        settings = self.set_dict

        Nx = settings["NX"]
        Nw = settings["NW"]
        L = settings["L"]
        W = settings["W"]
        Z_corner = settings["Z_CORNER"]

        def check_vectorise(x, name):
            """ Helper routine to perform sanity checks on `x`, followed by vectorisation """

            # Numeric types to check against
            scalar_types = (
                int, float, 
                np.int16, np.int32, np.int64, 
                np.float32, np.float64
            )

            # Loop over scalar types to check against type(x)
            scalar = False
            for t in scalar_types:
                scalar = scalar or isinstance(x, t)

            # If the input is a scalar: convert to vector
            if scalar:
                x = np.ones(Nw, dtype=float) * x

            # Check that the length of the vector equals Nw
            assert len(x) == Nw, f"The input vector `{name}` needs to be of length `NW` or be a scalar"
            
            return x

        dx = L / Nx  # Along-strike spacing
        if dw is None:
            dw = W / Nw  # Along-dip spacing

        # Perform sanity checks and create vectors (if needed)
        dip = check_vectorise(dip, name="dip")
        dw = check_vectorise(dw, name="dw")

        # Compute trigonometric function of the dip angle
        dip = np.deg2rad(dip)
        cd = np.cos(dip) * dw
        sd = np.sin(dip) * dw

        # Since the strike is constant, x always ranges from dx/2 to L-dx/2
        x = np.tile(np.linspace(0.5 * dx, L - 0.5 * dx, Nx), Nw)

        """
        The y and z vectors are the result of a summation, starting at the base
        Writing out the summation (starting at n = 0), we get:

        y[n] = y[0] + 0.5 * dw[0] * cos(dip[0]) + \sum_{k=1}^{n-1} dw[k] * cos(dip[k]) + 0.5 * dw[n] * cos(dip[n])
        z[n] = z[0] + 0.5 * dw[0] * sin(dip[0]) + \sum_{k=1}^{n-1} dw[k] * sin(dip[k]) + 0.5 * dw[n] * sin(dip[n])

        Instead of iterating, we use NumPy's cumsum function in a clever way to represent the partial sum
        """

        ccd = np.hstack([0, 0, np.cumsum(cd[1:-1])])
        csd = np.hstack([0, 0, np.cumsum(sd[1:-1])])

        y = cd[0] + 0.5 * cd + ccd
        y[0] = 0.5 * cd[0]

        z = sd[0] + 0.5 * sd + csd
        z[0] = 0.5 * sd[0]
        z += Z_corner

        def reshape_ravel(a, shape):
            return np.tile(a, shape).T.ravel()

        # Tiling (= repeating along-strike)
        shape = (Nx, 1)

        y = reshape_ravel(y, shape)
        z = reshape_ravel(z, shape)

        mesh_dict["X"] = x
        mesh_dict["Y"] = y
        mesh_dict["Z"] = z
        mesh_dict["DIP_W"] = reshape_ravel(np.rad2deg(dip), shape)
        mesh_dict["DW"] = reshape_ravel(dw, shape)

        pass

    def write_input(self):

        # The mesh should not have been rendered at this point
        assert self.mesh_rendered, "The mesh has not yet been rendered. Call render_mesh() before write_input()"

        # Optionally, output can be created to an external working directory
        # This is still experimental...
        if self.work_dir != "":
            if not os.path.isdir(self.work_dir):
                print("Creating %s" % (self.work_dir))
                subprocess.run(["mkdir", self.work_dir])
            print("Switching working directory to %s" % (self.work_dir))
            os.chdir(self.work_dir)

        settings = self.set_dict
        mesh = self.mesh_dict
        Nprocs = settings["NPROC"]
        delimiter = "    "

        # Check if multiple processors have been assigned for 0D and 1D faults
        assert not ((settings["MESHDIM"] != 2) and (settings["NPROC"] != 1)), 'Parallel processing (settings["NPROC"] > 1) is compatible only with settings["MESHDIM"] == 2'

        # Warn the user for legacy code
        if ((settings["VERBOSE"] == 1) and ("NTOUT" in settings.keys())):
            print("Warning: the setting `NTOUT` has been deprecated. Specify `NTOUT_OT`, `NTOUT_OX`, or `NTOUT_LOG` instead")

        # Define chunk size for each processor
        nwLocal = (settings["NW"] // Nprocs) * np.ones(Nprocs, dtype=int)
        nwLocal[Nprocs-1] += settings["NW"] % Nprocs
        nnLocal = 0

        # Number of elements along-strike
        Nx = settings["NX"]

        # Number of creep mechanisms
        N_creep = settings["SET_DICT_CNS"]["N_CREEP"]

        # Recalculate number of faults if necessary
        settings["N_FAULTS"] = len(np.unique(mesh["FAULT_LABEL"]))
        
        # If RSF: convert initial velocity to stress
        if (settings["FRICTION_MODEL"] == "RSF") and all(mesh["TAU"] < 0):
            if self.set_dict["VERBOSE"] == 1:
                print("Stress field not set, converting V0 into tau...")
            self.convert_velocity()

        # Loop over processor nodes
        for iproc in range(Nprocs):
            nloc = nwLocal[iproc]*settings["NX"]	# number of fault elements hosted on node
            iloc = np.zeros(nloc, dtype=int)		# buffer for indices of those elements

            if Nprocs == 1: filename = "qdyn.in"
            else: filename = "qdyn%06i.in" % (iproc)

            # Start building the contents of our QDYN input file (input_str)
            input_str = ""
            input_str += "%u%s meshdim\n" % (settings["MESHDIM"], delimiter)

            # Input specific to 3D faults
            if settings["MESHDIM"] == 2:
                input_str += "%u %u%s NX, NW\n" % (Nx, nwLocal[iproc], delimiter)
                input_str += "%.15g %.15g %.15g%s L, W, Z_CORNER\n" % (settings["L"], settings["W"], settings["Z_CORNER"], delimiter)
                for i in range(nwLocal[iproc]):
                    n = (i + np.sum(nwLocal[:iproc])) * Nx
                    input_str += "%.15g %.15g \n" % (mesh["DW"][n], mesh["DIP_W"][n])
            else:
                input_str += "%u%s NN\n" % (settings["N"], delimiter)
                input_str += "%.15g %.15g%s L, W\n" % (settings["L"], settings["W"], delimiter)

            # Input specific to 2D faults
            if settings["MESHDIM"] == 1:
                input_str += "%u%s finite\n" % (settings["FINITE"], delimiter)

            # Input specific to fault rheological model (CNS or RSF)
            if settings["FRICTION_MODEL"] == "RSF":
                input_str += "%u%s itheta_law\n" % (settings["SET_DICT_RSF"]["THETA_LAW"], delimiter)
                input_str += "%u%s i_rns_law\n" % (settings["SET_DICT_RSF"]["RNS_LAW"], delimiter)
            else: # CNS model
                input_str += "%u%s itheta_law\n" % (1, delimiter)
                input_str += "%u%s i_rns_law\n" % (3, delimiter)

            # Some more general settings

            # Add restart time
            input_str += "%.15g %u%s restart time step\n" % (settings["RESTART_TIME"], settings["RESTART_STEP"], delimiter)

            # Number of faults
            input_str += "%.15g%s fault number\n" % (settings["N_FAULTS"], delimiter)

            # If the CNS model is used, define the number of creep mechanisms
            if settings["FRICTION_MODEL"] == "CNS":
                # Raise an error when the default value has not been overwritten
                assert not (N_creep == -1), "The number of creep mechanisms was not set properly. This value should be determined in render_mesh()"
                input_str += "%u%s N_creep\n" % (N_creep, delimiter)

            input_str += "%u %u %u %u %u%s stress_coupling, thermal press.,fluid diffusion, variable k, localisation\n" % (settings["FEAT_STRESS_COUPL"], settings["FEAT_TP"],settings["FEAT_FLUID_DIFF"],settings["FEAT_VAR_K"], settings["FEAT_LOCALISATION"], delimiter)
            input_str += "%u %u %u %u %u %u %u %u %u %u%s ntout_log, ntout_ot, ntout_ox, nt_coord, nxout, nwout, nxout_DYN, nwout_DYN, ox_seq, ox_DYN\n" % (settings["NTOUT_LOG"], settings["NTOUT_OT"], settings["NTOUT_OX"], settings["IC"]+1, settings["NXOUT_OX"], settings["NWOUT_OX"], settings["NXOUT_DYN"], settings["NWOUT_DYN"], settings["OX_SEQ"], settings["OX_DYN"], delimiter)
            input_str += "%.15g %.15g %.15g %.15g %.15g %.15g%s beta, smu, lambda, v_th\n" % (settings["VS"], settings["MU"], settings["LAM"], settings["D"], settings["HD"], settings["V_TH"], delimiter)
            input_str += "%.15g %.15g%s Tper, Aper\n" % (settings["TPER"], settings["APER"], delimiter)
            input_str += "%.15g %.15g %.15g %.15g%s dt_try, dtmax, tmax, accuracy\n" % (settings["DTTRY"] ,settings["DTMAX"] ,settings["TMAX"] ,settings["ACC"] , delimiter)
            input_str += "%u%s nstop\n" % (settings["NSTOP"], delimiter)
            input_str += "%u %u%s DYN_FLAG, DYN_SKIP\n" % (settings["DYN_FLAG"], settings["DYN_SKIP"], delimiter)
            input_str += "%.15g %.15g %.15g%s M0, DYN_th_on, DYN_th_off\n" % (settings["DYN_M"], settings["DYN_TH_ON"], settings["DYN_TH_OFF"], delimiter)
            input_str += "%i %i%s FAULT_TYPE, SOLVER\n" % (settings["FAULT_TYPE"], settings["SOLVER"], delimiter)

            # Loop over all fault segments that are hosted on this processor node
            for i in range(nwLocal[iproc]):
                for j in range(settings["NX"]):
                    n = i*settings["NX"] + j
                    iloc[n] = j + i*settings["NX"] + nnLocal

            # Check if injection is required
            if settings["FEAT_FLUID_DIFF"] == 1:
                # Write the number of source 
                input_str += "%u%s nb_source\n" % (settings["SET_DICT_FLUID_DIFF"]["nb_source"], delimiter)
                
                # If the number of source is greater than one
                if settings["SET_DICT_FLUID_DIFF"]["nb_source"]>=1:
                    # Write the injection rate
                    for j in range(settings["SET_DICT_FLUID_DIFF"]["nb_source"]):
                        input_str += "%.15g %s" % (settings["SET_DICT_FLUID_DIFF"]["rate_injection"][j], delimiter)
                    # End injection 
                    input_str += "rate_injection %s\n" % (delimiter)
                    
                    # Write the index of injection
                    for j in range(settings["SET_DICT_FLUID_DIFF"]["nb_source"]):
                        input_str += "%.15g %s" % (settings["SET_DICT_FLUID_DIFF"]["index_injection"][j], delimiter)
                    # End injection 
                    input_str += "index_injection %s\n" % (delimiter)
                    
                    # Write the beginning time of injection
                    for j in range(settings["SET_DICT_FLUID_DIFF"]["nb_source"]):
                        input_str += "%.15g %s" % (settings["SET_DICT_FLUID_DIFF"]["t_injection_beg"][j], delimiter)
                    # End injection 
                    input_str += "t_injection_beg %s\n" % (delimiter)
                    
                    # Write the end time of injection
                    for j in range(settings["SET_DICT_FLUID_DIFF"]["nb_source"]):
                        input_str += "%.15g %s" % (settings["SET_DICT_FLUID_DIFF"]["t_injection_end"][j], delimiter)
                    # End injection 
                    input_str += "t_injection_end %s\n" % (delimiter)

                
                
                



            # Check for fault rheology model
            if settings["FRICTION_MODEL"] == "CNS":
                for i in range(nloc):
                    # Basic parameters (initial conditions, granular flow, porosity)
                    input_str += "%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g " % (mesh["SIGMA"][iloc[i]], mesh["TAU"][iloc[i]], mesh["PHI_INI"][iloc[i]], mesh["V_PL"][iloc[i]], mesh["A_TILDE"][iloc[i]], mesh["MU_TILDE_STAR"][iloc[i]], mesh["Y_GR_STAR"][iloc[i]], mesh["H"][iloc[i]], mesh["PHI_C"][iloc[i]], mesh["PHI0"][iloc[i]])
                    # Creep parameters (3 for each creep mechanism)
                    for j in range(N_creep):
                        input_str += "%.15g %.15g %.15g " % (mesh["A%i" % j][iloc[i]], mesh["N%i" % j][iloc[i]], mesh["M%i" % j][iloc[i]])
                    # Last couple of parameters
                    input_str += "%.15g %.15g %.15g\n" % (mesh["THICKNESS"][iloc[i]], mesh["IOT"][iloc[i]], mesh["IASP"][iloc[i]])
            else: # RSF model
                for i in range(nloc):
                    input_str += "%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n" % (mesh["SIGMA"][iloc[i]], mesh["TAU"][iloc[i]], mesh["TH_0"][iloc[i]], mesh["A"][iloc[i]], mesh["B"][iloc[i]], mesh["DC"][iloc[i]], mesh["V1"][iloc[i]], mesh["V2"][iloc[i]], mesh["MU_SS"][iloc[i]], mesh["V_SS"][iloc[i]], mesh["IOT"][iloc[i]], mesh["IASP"][iloc[i]], mesh["CO"][iloc[i]], mesh["V_PL"][iloc[i]], mesh["INV_VISC"][iloc[i]])

            # Check if localisation is requested
            print('--------')
            print(settings["FEAT_LOCALISATION"])
            print(settings["FRICTION_MODEL"])
            print('--------')
            
            
            
            if settings["FEAT_LOCALISATION"] == 1:
                if settings["FRICTION_MODEL"] != "CNS":
                    raise ValueError("Localisation is compatible only with the CNS friction model")
                for i in range(nloc):
                    # Basic parameters (degree of localisation, initial bulk porosity)
                    input_str += "%.15g %.15g " % (mesh["LOCALISATION"][iloc[i]], mesh["PHI_INI_BULK"][iloc[i]])
                    # Creep parameters (3 for each creep mechanism)
                    for j in range(N_creep):
                        input_str += "%.15g %.15g %.15g " % (mesh["A%i_bulk" % j][iloc[i]], mesh["N%i_bulk" % j][iloc[i]], mesh["M%i_bulk" % j][iloc[i]])
                    # End the line
                    input_str += "\n"

            # Check if thermal pressurisation is requested
            if settings["FEAT_TP"] == 1:
                for i in range(nloc):
                    input_str += "%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n" % (mesh["RHOC"][iloc[i]], mesh["BETA"][iloc[i]], mesh["ETA"][iloc[i]], mesh["HALFW"][iloc[i]], mesh["K_T"][iloc[i]], mesh["K_P"][iloc[i]], mesh["LAM"][iloc[i]], mesh["P_A"][iloc[i]], mesh["T_A"][iloc[i]], mesh["DILAT_FACTOR"][iloc[i]])
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             # Check if fluid diffusion is requested
            if settings["FEAT_FLUID_DIFF"] == 1:
                 for i in range(nloc):
                     input_str += "%.15g %.15g %.15g %.15g %.15g %.15g\n" % (mesh["RHOF"][iloc[i]], mesh["BETA"][iloc[i]], mesh["ETA"][iloc[i]], mesh["PHI"][iloc[i]], mesh["PERMEABILITY"][iloc[i]], mesh["P_A"][iloc[i]])

            # Check if variable permeability is requested
            if settings["FEAT_VAR_K"] == 1:
                for i in range(nloc):
                    input_str += "%.15g %.15g %.15g %.15g %.15g\n" % (mesh["Tk"][iloc[i]], mesh["Lk"][iloc[i]], mesh["Snk"][iloc[i]], mesh["kmin"][iloc[i]], mesh["kmax"][iloc[i]])



            
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            # Add mesh grid location information + restart_slip values
            for i in range(nloc):
                input_str += "%.15g %.15g %.15g %.15g %.15g %.15g\n" % (mesh["X"][iloc[i]], mesh["Y"][iloc[i]], mesh["Z"][iloc[i]], mesh["DIP_W"][iloc[i]], mesh["FAULT_LABEL"][iloc[i]], mesh["RESTART_SLIP"][iloc[i]])

            nnLocal += settings["NX"]*nwLocal[iproc]

            # Write QDYN input file
            with open(filename, "w") as f:
                f.write(input_str)

        self.qdyn_input_written = True

        return True

    # Run QDYN
    def run(self, run=True, test=False, unit=False, restart=False):

        assert self.qdyn_input_written or unit, "Input file has not yet been written. First call write_input() to render QDyn input file"

        # Executable file (including path)
        qdyn_exec = os.path.join(self.qdyn_path, "qdyn")

        # If we're on a Windows 10 + Unix subsystem, we need to adjust executable path
        if self.W10_bash is True:
            # Replace Win-style backward slashes with Unix forward slashes
            qdyn_exec_bash = qdyn_exec.replace("\\", "/")
            # Find drive letter (usually C:, but not always)
            drive = re.match(r"([A-Z](?=:))", qdyn_exec_bash).group()
            # Replace drive with mounted partition
            qdyn_exec_bash = qdyn_exec_bash.replace("%s:" % (drive), "/mnt/%s" % (drive.lower()))

        if restart:

            if self.__dict__.get("ox_last") is None:
                self.read_output(read_ot=False, read_ox=False, read_ox_dyn=False, read_ox_last=True)
            
            assert self.mesh_rendered, "The mesh has not yet been rendered"

            self.mesh_dict["TAU"] = self.ox_last["tau"].values
            self.mesh_dict["SIGMA"] = self.ox_last["sigma"].values
            if self.set_dict["FRICTION_MODEL"] == "RSF":
                self.mesh_dict["TH_0"] = self.ox_last["theta"].values
            elif self.set_dict["FRICTION_MODEL"] == "CNS":
                self.mesh_dict["PHI_INI"] = self.ox_last["theta"].values
            
            if self.set_dict["FEAT_TP"] == 1:
                self.mesh_dict["P_A"] = self.ox_last["P"].values
                self.mesh_dict["T_A"] = self.ox_last["T"].values

            self.mesh_dict["RESTART_SLIP"] = self.ox_last["slip"].values
            self.set_dict["RESTART_TIME"] = self.ox_last["t"].iloc[0]
            self.set_dict["RESTART_STEP"] = self.ox_last["step"].iloc[0]

            if self.set_dict["VERBOSE"] == 1:
                print("Rewriting input file for simulation restart")

            self.write_input()

        flags = []

        if test is True:
            self.set_dict["VERBOSE"] = 0

        if self.set_dict["DEBUG"] == 1:
            flags.append("debug")
        if self.set_dict["VERBOSE"] == 1:
            flags.append("verbose")
        if restart:
            flags.append("restart")
        if unit is True:
            flags.append("test")

        # If serial (NPROC = 1)
        if self.set_dict["NPROC"] == 1:

            # If we're on a Windows 10 + Unix subsystem, call bash
            if self.W10_bash is True: cmd = ["bash", "-c", "\"%s\"" % (qdyn_exec_bash)]
            # If we're on Unix, simply call the qdyn executable directly
            else: cmd = [qdyn_exec]
            
        else: # MPI parallel

            MPI_path = self.set_dict["MPI_PATH"]

            # If we're on a Windows 10 + Unix subsystem, call bash
            if self.W10_bash is True:
                cmd = ["bash", "-c",
                "\"%s -np %i %s\"" % (MPI_path, self.set_dict["NPROC"], qdyn_exec_bash)]
            # If we're on Unix, simply call the mpi executable directly
            else:
                cmd = [MPI_path, "-np", "%i" % (self.set_dict["NPROC"]), qdyn_exec]

        # Add command line flags
        if len(flags) > 0:
            cmd += flags

        
        # Sometimes it is convenient to launch run(restart=True) without calling
        # the executable (for instance on a cluster)
        if not run:
            return None
        
        # Run QDYN
        result = subprocess.run(cmd)
        if unit: return result.returncode

        # If a suffix is requested, rename output files
        suffix = self.set_dict["SUFFIX"]
        if suffix != "":
            for filename in os.listdir("."):
                if filename.startswith("fort") and "_" not in filename: os.rename(filename, "%s_%s" % (filename, suffix))

        return True


    # Read QDYN output data
    # This modified function allows to retrieve the outputs from a folder other than the one of the python script
    # It can also read the output of the last snapshot

    def read_output(self, mirror=False, path_output=None, read_ot=True, filename_ot="output_ot",
                    filename_vmax="output_vmax", read_ox=True,
                    filename_ox="output_ox", read_ox_dyn=False, filename_ox_last="output_ox_last", read_ox_last=False,
                    filename_fault="output_fault", read_fault=False):

        # Output file contents depends on the requested features
        ## Time series (ot)
        nheaders_ot = 5
        quants_ot = ("step", "t", "potcy", "pot_rate", "v", "theta", "tau", "dtau_dt", "slip", "sigma")
        if self.set_dict["FEAT_TP"] == 1:
            nheaders_ot += 1
            quants_ot += ("P", "T")
        if self.set_dict["FEAT_FLUID_DIFF"] == 1:
             nheaders_ot += 1
             quants_ot += ("P",)
        quants_ot += ("fault_label",)

        # Vmax
        nheaders_vmax = 3
        quants_vmax = ("step", "t", "ivmax", "v", "theta", "tau", "dtau_dt", "slip", "sigma")
        if self.set_dict["FEAT_TP"] == 1:
            nheaders_vmax += 1
            quants_vmax += ("P", "T")
            
        if self.set_dict["FEAT_FLUID_DIFF"] == 1:
            nheaders_vmax += 1
            quants_vmax += ("P",)
        
            
        quants_vmax += ("fault_label",)

        # Standard snapshots (ox)
        quants_ox = ("step", "t", "x", "y", "z", "v", "theta", "tau", "tau_dot", "slip", "sigma")
        if self.set_dict["FEAT_TP"] == 1:
            quants_ox += ("P", "T")
            
        if self.set_dict["FEAT_FLUID_DIFF"] == 1:
            quants_ox += ("P",)
            
        quants_ox += ("fault_label",)
        
        # Fault time-series
        nheaders_fault = 2
        quants_fault = ("step", "t", "potcy_fault", "pot_rate_fault", "vmax_fault", "ivmax_fault")

        # If time series data is requested
        if read_ot:
            iot = self.mesh_dict["IOT"]
            inds = (iot == 1)
            N_iot = inds.sum()
            self.N_iot = N_iot

            iot = np.arange(len(self.mesh_dict["IOT"]))[inds]
            self.iot_inds = iot
            self.ot = [None] * N_iot

            for n, i in enumerate(iot):

                # Check output directory
                if path_output!=None:
                    filename_ot = path_output + filename_ot

                filename_iot = "%s_%i" % (filename_ot, i)

                df = read_csv(
                    filename_iot, header=None, skiprows=nheaders_ot,
                    names=quants_ot, sep="\s+"
                )

                # Sanitise output (check for near-infinite numbers, etc.)
                df = df.apply(pd.to_numeric, errors="coerce")

                # Discard duplicate rows from duplicate time-steps
                df = df.drop_duplicates(subset=["step"], keep="first") 
                self.ot[n] = df

            # Check output directory
            if path_output!=None:
                filename_vmax = path_output + filename_vmax

            self.ot_vmax = read_csv(
                filename_vmax, header=None, skiprows=nheaders_vmax,
                names=quants_vmax, sep="\s+"
            )
            # Discard duplicate rows from duplicate time-steps
            self.ot_vmax = self.ot_vmax.drop_duplicates(subset=["step"], keep="first") 

        else:
            self.ot = None
            self.ot_vmax = None
            self.N_iot = 0

        ## Standard snapshots
        if read_ox:

            # Read snapshot output
            # Check output directory
            if path_output!= None:
                filename_ox= path_output + filename_ox

            data_ox = read_csv(filename_ox, header=None, names=quants_ox, sep="\s+", comment="#",low_memory=False)

            # Store snapshot data in self.ox
            self.ox = data_ox

            # Sanitise output (check for near-infinite numbers, etc.)
            self.ox = self.ox.apply(pd.to_numeric, errors="coerce")

            # Discard duplicate rows from duplicate time-steps
            self.ox = self.ox.drop_duplicates(subset=["t", "x", "y", "z"], keep="first") 

            # If free surface was generated manually (i.e. without FINITE = 2 or 3),
            # take only half data set (symmetric around first element)
            if mirror == True:
                data_ox = data_ox.loc[(data_ox["x"] > 0)]
                self.ox = data_ox
        else:
            self.ox = None

        ## Last snapshot
        if read_ox_last == True:

            # Read snapshot output
            # Check output directory
            if path_output!= None:
                filename_ox_last = path_output + filename_ox_last

            data_ox_last = read_csv(filename_ox_last, header=None, names=quants_ox, sep="\s+", comment="#")

            # Store snapshot data in self.ox
            self.ox_last = data_ox_last

            # Sanitise output (check for near-infinite numbers, etc.)
            self.ox_last = self.ox_last.apply(pd.to_numeric, errors="coerce")

            # If free surface was generated manually (i.e. without FINITE = 2 or 3),
            # take only half data set (symmetric around first element)
            if mirror == True:
                data_ox_last = data_ox_last.loc[(data_ox_last["x"] > 0)]
                self.ox_last = data_ox_last
        else:
            self.ox_last = None

        ## Dynamic snapshot
        if read_ox_dyn == True:

            # Check directory
            if path_output != None:
                ox_dyn_files_pre = np.array([file for file in os.listdir(path_output) if file.startswith("output_dyn_pre")])
                ox_dyn_files_post = np.array([file for file in os.listdir(path_output) if file.startswith("output_dyn_post")])
                ox_dyn_files_rup = np.array([file for file in os.listdir(path_output) if file.startswith("output_dyn_max")])
            else:                
                ox_dyn_files_pre = np.array([file for file in os.listdir(".") if file.startswith("output_dyn_pre")])
                ox_dyn_files_post = np.array([file for file in os.listdir(".") if file.startswith("output_dyn_post")])
                ox_dyn_files_rup = np.array([file for file in os.listdir(".") if file.startswith("output_dyn_max")])

            # Read snapshot output
            data_ox_dyn_pre = [None] * len(ox_dyn_files_pre)
            data_ox_dyn_post = [None] * len(ox_dyn_files_post)
            data_ox_dyn_rup = [None] * len(ox_dyn_files_rup)
            for i in range(len(ox_dyn_files_pre)):

                quants_rup = ("x", "y", "z", "t", "tau_max", "t_v_max", "v_max")

                # Read pre-rupture file
                data_ox_dyn_pre[i] = read_csv(
                    ox_dyn_files_pre[i], header=None, names=quants_ox,
                    sep="\s+", comment="#"
                )
                # Sanitise output (check for near-infinite numbers, etc.)
                data_ox_dyn_pre[i] = data_ox_dyn_pre[i].apply(pd.to_numeric, errors="coerce")
            
                # Discard duplicate rows from duplicate time-steps
                data_ox_dyn_pre[i] = data_ox_dyn_pre[i].drop_duplicates(subset=["t", "x", "y", "z"], keep="first")

                # Read pre-rupture file
                data_ox_dyn_post[i] = read_csv(
                    ox_dyn_files_post[i], header=None, names=quants_ox,
                    sep="\s+", comment="#"
                )
                # Sanitise output (check for near-infinite numbers, etc.)
                data_ox_dyn_post[i] = data_ox_dyn_post[i].apply(pd.to_numeric, errors="coerce")

                # Discard duplicate rows from duplicate time-steps
                data_ox_dyn_post[i] = data_ox_dyn_post[i].drop_duplicates(subset=["t", "x", "y", "z"], keep="first")

                # Rupture stats file
                data_ox_dyn_rup[i] = read_csv(
                    ox_dyn_files_rup[i], header=None, names=quants_rup,
                    sep="\s+", comment="#"
                )
                # Sanitise output (check for near-infinite numbers, etc.)
                data_ox_dyn_rup[i] = data_ox_dyn_rup[i].apply(pd.to_numeric, errors="coerce")

                # Discard duplicate rows from duplicate time-steps
                data_ox_dyn_rup[i] = data_ox_dyn_rup[i].drop_duplicates(subset=["t", "x", "y", "z"], keep="first") 

            # Store snapshot data in self.ox
            self.ox_dyn_pre = data_ox_dyn_pre
            self.ox_dyn_post = data_ox_dyn_post
            self.ox_dyn_rup = data_ox_dyn_rup


        # If fault data is requested

        if read_fault:
            
            # Number of faults
            N_faults = self.set_dict["N_FAULTS"]

            # List where each element is a DataFrame with the output of one fault
            self.fault = [None] * N_faults

            # Check output directory
            if path_output!=None:
                filename_fault = path_output + filename_fault 

            # Number of columns with fault output (excluding step/time)
            stride = 4
            N_cols = stride * N_faults

            # Counter of fault number for loop
            n = 0
            
            # Read the output file in groups of columns corresponding to each fault
            for i in range(2, N_cols + 2, stride):
                # Column list for one fault (time, potency, potency rate, vmax_fault and ivmax_fault)
                col_list = [0, 1] + list(range(i, i + stride))
            
                # Read output file
                self.fault[n] = read_csv(
                filename_fault, header=None, skiprows=nheaders_fault, usecols=col_list,
                names=quants_fault, sep="\s+"
                )

                # Discard duplicate rows from duplicate time-steps
                self.fault[n] = self.fault[n].drop_duplicates(subset=["step"], keep="first")
                n =+1
        return True


    # For compatibility with the viscosity model, the initial state of the
    # fault is expressed in terms of the shear stress (and not velocity).
    # For rate-and-state friction, the initial velocity needs to be converted
    # initial stress
    def convert_velocity(self):

        set_dict = self.set_dict["SET_DICT_RSF"]
        mesh_dict = self.mesh_dict

        # Step 1: load V0, theta0
        v = mesh_dict["V_0"]
        theta = mesh_dict["TH_0"]

        # Step 2: select friction law
        if set_dict["RNS_LAW"] == 0:
            fric_law = self.friction_regular
            v_law = self.v_regular
        elif set_dict["RNS_LAW"] == 1:
            fric_law = self.friction_cutoff
            v_law = self.v_cutoff
        elif set_dict["RNS_LAW"] == 2:
            fric_law = self.friction_regularised
            v_law = self.v_regularised
        else:
            print("Friction law %d not recognised" % set_dict["RNS_LAW"])
            exit()

        # Step 3: if viscosity == 0: use friction law directly
        if all(mesh_dict["INV_VISC"] == 0):
            tau = fric_law(v, theta) * (mesh_dict["SIGMA"] - mesh_dict["P_A"])

        # Step 4: if viscosity > 0: solve for tau -> v0 = vfric + vvisc
        else:
            tau = self.solve_friction(v, theta, v_law)

        # Step 5: store result in mesh_dict
        self.mesh_dict["TAU"] = tau

        return True

    def friction_regular(self, v, theta):
        mesh_dict = self.mesh_dict
        mu_star = mesh_dict["MU_SS"]
        v_star = mesh_dict["V_SS"]
        a = mesh_dict["A"]
        b = mesh_dict["B"]
        Dc = mesh_dict["DC"]
        return mu_star + a * np.log(v / v_star) + b * np.log(theta * v_star / Dc)

    def v_regular(self, tau, theta):
        mesh_dict = self.mesh_dict
        mu = tau / (mesh_dict["SIGMA"] - mesh_dict["P_A"])
        mu_star = mesh_dict["MU_SS"]
        v_star = mesh_dict["V_SS"]
        a = mesh_dict["A"]
        b = mesh_dict["B"]
        Dc = mesh_dict["DC"]
        return v_star * np.exp((mu - mu_star - b * np.log(theta * v_star / Dc)) / a)

    def friction_cutoff(self, v, theta):
        mesh_dict = self.mesh_dict
        mu_star = mesh_dict["V_SS"]
        v1 = mesh_dict["V1"]
        v2 = mesh_dict["V2"]
        a = mesh_dict["A"]
        b = mesh_dict["B"]
        Dc = mesh_dict["DC"]
        return mu_star - a * np.log(v1 / v + 1) + b * np.log(theta * v2 / Dc + 1)

    def v_cutoff(self, tau, theta):
        mesh_dict = self.mesh_dict
        mu = tau / (mesh_dict["SIGMA"] - mesh_dict["P_A"])
        mu_star = mesh_dict["V_SS"]
        v1 = mesh_dict["V1"]
        v2 = mesh_dict["V2"]
        a = mesh_dict["A"]
        b = mesh_dict["B"]
        Dc = mesh_dict["DC"]
        return v1 / (np.exp(-(mu - mu_star - b * np.log(theta * v2 / Dc + 1)) / a) - 1)

    def friction_regularised(self, v, theta):
        mesh_dict = self.mesh_dict
        mu_star = mesh_dict["MU_SS"]
        v_star = mesh_dict["V_SS"]
        a = mesh_dict["A"]
        b = mesh_dict["B"]
        Dc = mesh_dict["DC"]

        exp_term = np.exp((mu_star + b * np.log(theta * v_star / Dc)) / a)

        return a * np.arcsinh(v * exp_term / (2 * v_star))

    def v_regularised(self, tau, theta):
        mesh_dict = self.mesh_dict
        mu = tau / (mesh_dict["SIGMA"] - mesh_dict["P_A"])
        mu_star = mesh_dict["MU_SS"]
        v_star = mesh_dict["V_SS"]
        a = mesh_dict["A"]
        b = mesh_dict["B"]
        Dc = mesh_dict["DC"]

        exp_term = np.exp(-(mu_star + b * np.log(theta * v_star / Dc)) / a)

        return 2 * v_star * np.sinh(mu / a) * exp_term

    def solve_friction(self, v, theta, v_law):
        inv_visc = self.mesh_dict["INV_VISC"]
        x0 = 0.5 * (self.mesh_dict["SIGMA"] - self.mesh_dict["P_A"])

        root_func = lambda x: v - v_law(x, theta) - x * inv_visc
        result = root(root_func, x0, method="lm")
        tau = result["x"]

        assert all(np.isfinite(tau)), "The friction could not be estimated"

        return tau

