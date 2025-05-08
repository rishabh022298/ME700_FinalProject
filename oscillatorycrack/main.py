from dolfin import *
from mshr import *
import numpy as np
import time

start_time = time.time()

# Domain and mesh
L = 10
h = 2

domain = Rectangle(Point(-L, -h/2), Point(L, h/2))
domain.set_subdomain(1, Rectangle(Point(-L,-0.01), Point(0.5,0.0)))  # For crack
mesh = generate_mesh(domain, 100)
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim(), mesh.domains())
boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)
dx = Measure('dx', domain=mesh, subdomain_data=subdomains)

File("./ResultsDir/subdomains.pvd") << subdomains

def Crack(x):
    return abs(x[1]) < 1e-02 and x[0] <= 0.5

r = Expression(('x[0]','x[1]'),degree=1)

# Function Spaces
V = FunctionSpace(mesh, 'CG', 1)
W = VectorFunctionSpace(mesh, 'CG', 1)
WW = FunctionSpace(mesh, 'DG', 0)
p, q = TrialFunction(V), TestFunction(V)
FS = FunctionSpace(mesh, 'Lagrange', 1)
P1 = FiniteElement('Lagrange', mesh.ufl_cell(), degree=2)
R1 = FiniteElement('Real', mesh.ufl_cell(), 0)
MFS = FunctionSpace(mesh, MixedElement([(P1*P1),(R1*R1),R1]))
f = Function(MFS)
v_f = TestFunction(MFS)
u, c_trans, c_rot = split(f)

# Material and system parameters
## Change values here to check different regions in the phase field
Gc =  1          # Critical energy release rate
l = 0.015
lmbda = 1000000
mu = 800000
P = Constant(5.0)  
Diffusion_const = Constant(1.0)
vel = Constant(P*Diffusion_const/(h/2))
alpha = 0.01        # Coefficient of thermal expansion

# Moving temperature field
Temp_expr = Expression("x[0] >= vel*t ? (1.0 - exp(-P * (x[0] - vel*t))) : 0.0", degree=2, P=P, t=0, vel = vel)
# Constitutive functions
def epsilon(u):
    return sym(grad(u)) - alpha*as_matrix([[Temp_expr,0],[0,Temp_expr]])
def sigma(u):
    return (2.0*mu*epsilon(u) + lmbda*tr(epsilon(u))*Identity(len(u)))
def psi(u):
    return (0.5 * (lmbda + mu) * (0.5 * (tr(epsilon(u)) + abs(tr(epsilon(u)))))**2 + mu * inner(dev(epsilon(u)), dev(epsilon(u))))		
def H(uold, unew, Hold):
    return conditional(lt(psi(uold), psi(unew)), psi(unew), Hold)
def bulk_energy_density(u):
    return 0.5*inner(sigma(u),epsilon(u))

# Boundary conditions
bc_phi = [DirichletBC(V, Constant(1.0), Crack)]
bc_u = []

# Variational form
unew, uold = Function(W), Function(W)
pnew, pold, Hold = Function(V), Function(V), Function(V)

# Bounds on damage parameter
lb = Function(V)
ub = Function(V)
lb.vector()[:] = 0.0  # Lower bound: p ≥ 0
ub.vector()[:] = 1.0  # Upper bound: p ≤ 1

# Energy functionals
E_du = ((1.0 - pold)**2) * bulk_energy_density(u) * dx + dot(c_trans,u)*dx + c_rot*(r[0]*u[1]-r[1]*u[0])*dx
E_phi = (((l**2) * dot(grad(p), grad(q))) + ((2*l/Gc) * H(unew, uold, Hold) +1 ) * p * q )* dx - (2*l/Gc) * H(unew, uold, Hold) * q * dx
E_phi = action(E_phi, pnew)

# Residuals and setting up the solver
Res_u = derivative(E_du, f, v_f)
Jac_u = derivative(Res_u, f) 
p_disp = NonlinearVariationalProblem(Res_u, f, bc_u, Jac_u)
solver_disp = NonlinearVariationalSolver(p_disp)

Jac_p = derivative(E_phi,pnew,p)
p_phi = NonlinearVariationalProblem(E_phi, pnew, bc_phi, Jac_p)

p_phi.set_bounds(lb, ub)

solver_phi = NonlinearVariationalSolver(p_phi)

# Solver parameters
snes_Rtol         = 1e-6 # relative tolerance for SNES solver (phase field eq.)
snes_Atol         = 1e-6 # absolute tolerance for SNES solver (phase field eq.)
snes_maxiter      = 300  # max. iteration for SNEs solver (phase field eq.)

snes_prm = {"nonlinear_solver": "snes",
                "snes_solver"     : { "method": "vinewtonssls",
                                    "line_search": "basic",
                                    "maximum_iterations": snes_maxiter,
                                    "relative_tolerance": snes_Rtol,
                                    "absolute_tolerance": snes_Atol,
                                    "report": True,
                                    "error_on_nonconvergence": False,
                                    }}
solver_phi.parameters.update(snes_prm)

# Initialization and simulation loop
t = 0

deltaT = 0.001  # Load steps
tol = 4e-3      # For staggered scheme
conc_f = File("./ResultsDir/phi.pvd")
conc_f << pnew
conc_u = File("./ResultsDir/u.pvd")
conc_u << unew

# Visualizing moving temperature field
Temp_space = FunctionSpace(mesh, 'CG', 1)
Temp_func = Function(Temp_space)
temp_file = File("./ResultsDir/temperature.pvd")

# Stresses and strains
fileSxx = File("./ResultsDir/sigma_xx.pvd")
fileSyy = File("./ResultsDir/sigma_yy.pvd")
fileSxy = File("./ResultsDir/sigma_xy.pvd")

fileexx = File("./ResultsDir/epsilon_xx.pvd")
fileeyy = File("./ResultsDir/epsilon_yy.pvd")
fileexy = File("./ResultsDir/epsilon_xy.pvd")

while t <= 0.1:   
    if t<0.020:
        t += deltaT
    else:
        t += deltaT/2          # Smaller load steps for stable crack propagation
    Temp_expr.t = t
    iter, err = 0, 1

    while err > tol:
        iter += 1
        
        solve(Res_u == 0, f, bc_u,
              solver_parameters={
                  "newton_solver": {
                      "maximum_iterations": 100,
                      "relative_tolerance": 1e-6,
                      "absolute_tolerance": 1e-8,
                      "linear_solver": "mumps",  # or "lu"
                      "report": True,
                      "error_on_nonconvergence": True,
                      "relaxation_parameter": 0.4
                  }
              })
        
        u_mixed, _, _ = f.split(deepcopy=True)
        unew.assign(u_mixed)
        solver_phi.solve()
        err_u = errornorm(unew, uold, norm_type='l2')
        err_phi = errornorm(pnew, pold, norm_type='l2')
        print(f"t = {t:.5f}, iter = {iter}, err_u = {err_u:.3e}, err_phi = {err_phi:.3e}")
        err = max(err_u, err_phi)

        uold.assign(unew)
        pold.assign(pnew)
        Hold.assign(project(psi(unew), WW))

        if err < tol:
            print('Iterations:', iter, ', Total time', t)
            conc_f << pnew
            conc_u << unew
            
            sigma_xx = project(sigma(unew)[0,0], FS)
            sigma_xx.rename("sigma_xx","")
            fileSxx << sigma_xx
            sigma_yy = project(sigma(unew)[1,1], FS)
            sigma_yy.rename("sigma_yy","")
            fileSyy << sigma_yy
            sigma_xy = project(sigma(unew)[0,1], FS)
            sigma_xy.rename("sigma_xy","")
            fileSxy << sigma_xy
            epsilon_xx = project(epsilon(unew)[0,0], FS)
            epsilon_xx.rename("epsilon_xx","")
            fileexx << epsilon_xx
            epsilon_yy = project(epsilon(unew)[1,1], FS)
            epsilon_yy.rename("epsilon_yy","")
            fileeyy << epsilon_yy
            epsilon_xy = project(epsilon(unew)[0,1], FS)
            epsilon_xy.rename("epsilon_xy","")
            fileexy << epsilon_xy

            Temp_func.assign(project(Temp_expr, Temp_space))
            temp_file << (Temp_func, t)

            # Generating text files for each load step to mark the location of points in the domain with phase field parameter value > 0.95
            # User may comment the above file generation lines (or below) to save some time.
            dof_coords = V.tabulate_dof_coordinates().reshape((-1, 2))
            p_vals = pnew.vector().get_local()
            mask = p_vals > 0.95
            np.savetxt(f"./ResultsDir/phasefield_above_095_t{t:.4f}.txt", dof_coords[mask], fmt="%.6f", header="x y")
            
end_time = time.time()

execution_time = end_time - start_time
print("Time for completion:", execution_time, "\n")
print('Simulation completed')
