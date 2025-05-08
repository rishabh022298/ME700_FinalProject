from dolfin import *
from mshr import *
import numpy as np

# Define radius of regions to calculate J Integral
R1 = 0.25
R2 = 0.3

# Dimensions of the domain
L = 1
H = 1

# Mesh
domain = Rectangle(Point(-L/2, -H/2), Point(L/2, H/2))
domain.set_subdomain(1, (Circle(Point(0, 0), R2) - Circle(Point(0, 0), R1) - Rectangle(Point(-R2-0.02, -0.01), Point(-R1+0.02,0.01)))) # domain for J-integral calculations
domain.set_subdomain(2, Rectangle(Point(-L/2,-0.01), Point(0,0)))    # for setting up the crack
mesh = generate_mesh(domain, 200)
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim(), mesh.domains())
boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)
dx = Measure('dx', domain=mesh, subdomain_data=subdomains)

# Save markers for visualization
File("subdomains.pvd") << subdomains

# Define boundary classes
class TopBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.5)

class BottomBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], -0.5)

def Crack(x):
    return abs(x[1]) < 1e-03 and x[0] <= 0.0

# Mark boundaries
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundaries.set_all(0)

top_boundary = TopBoundary()
bottom_boundary = BottomBoundary()

top_boundary.mark(boundaries, 3)      
bottom_boundary.mark(boundaries, 4)   

ds = Measure("ds", subdomain_data=boundaries)
n = FacetNormal(mesh)

# Function Spaces
V = FunctionSpace(mesh, 'CG', 1)
W = VectorFunctionSpace(mesh, 'CG', 1)
WW = FunctionSpace(mesh, 'DG', 0)
p, q = TrialFunction(V), TestFunction(V)
u, v = TrialFunction(W), TestFunction(W)
FS = FunctionSpace(mesh, 'Lagrange', 1)

# Material parameters
Gc =  2.7
l = 0.015
lmbda = 121.1538e3
mu = 80.7692e3

# Constitutive functions
def epsilon(u):
    return sym(grad(u))
def sigma(u):
    return 2.0*mu*epsilon(u) + lmbda*tr(epsilon(u))*Identity(len(u))
def psi(u):
    return 0.5 * (lmbda + mu) * (0.5 * (tr(epsilon(u)) + abs(tr(epsilon(u)))))**2 + mu * inner(dev(epsilon(u)), dev(epsilon(u)))		
def H(uold, unew, Hold):
    return conditional(lt(psi(uold), psi(unew)), psi(unew), Hold)

# Define boundary conditions
load_top = Expression("t", t=0.0, degree=1)
load_bot = Expression("t", t=0.0, degree=1)

bcbot = DirichletBC(W.sub(1), load_bot, bottom_boundary)
bctop = DirichletBC(W.sub(1), load_top, top_boundary)
bc_u = [bcbot, bctop]

bc_phi = [DirichletBC(V, Constant(1.0), Crack)]

# Variational form
unew, uold = Function(W), Function(W)
pnew, pold, Hold = Function(V), Function(V), Function(V)

E_du = ((1.0 - pold)**2) * inner(grad(v), sigma(u)) * dx
E_phi = (Gc * l * inner(grad(p), grad(q)) + ((Gc / l) + 2.0 * H(uold, unew, Hold)) * inner(p, q) - 2.0 * H(uold, unew, Hold) * q) * dx

p_disp = LinearVariationalProblem(lhs(E_du), rhs(E_du), unew, bc_u)
p_phi = LinearVariationalProblem(lhs(E_phi), rhs(E_phi), pnew, bc_phi)

solver_disp = LinearVariationalSolver(p_disp)
solver_phi = LinearVariationalSolver(p_phi)

# Weight function for domain J integral
q_func_space = FunctionSpace(mesh, 'CG', 1)
q_expr = Expression('sqrt(pow(x[0], 2) + pow(x[1], 2)) < r1 ? \
                     1.0 : (sqrt(pow(x[0], 2) + pow(x[1], 2)) > r2 ? 0.0 : \
                     (r2 - sqrt(pow(x[0], 2) + pow(x[1], 2))) / (r2 - r1))', \
                          degree=1, r1=R1, r2=R2)
q_J_integral = interpolate(q_expr, q_func_space)

# Initialization and simulation loop
t = 0
u_r = 0.007
deltaT = 0.1
tol = 1e-3
conc_f = File("./ResultsDir/phi.pvd")
fname2 = open('JInt.txt', 'w')
se_file = File("./ResultsDir/se.pvd")    # Strain Energy
while t <= 0.3:
    t += deltaT
    load_top.t = t * u_r
    load_bot.t = -t * u_r
    iter, err = 0, 1

    while err > tol:
        iter += 1
        solver_disp.solve()
        solver_phi.solve()
        err_u = errornorm(unew, uold, norm_type='l2')
        err_phi = errornorm(pnew, pold, norm_type='l2')
        err = max(err_u, err_phi)

        uold.assign(unew)
        pold.assign(pnew)
        Hold.assign(project(psi(unew), WW))

        if err < tol:
            print('Iterations:', iter, ', Total time', t)

            if round(t * 1e4) % 10 == 0:
                conc_f << pnew
                se_func = project(0.5*inner(epsilon(unew),sigma(unew)),V)
                se_func.rename("se_func","strain energy")
                se_file << se_func

                # Compute J-integral
                J_integral = assemble(dot((sigma(unew)*(grad(unew)) - 0.5*inner(sigma(unew),epsilon(unew))*Identity(2))[:,0] , grad(q_J_integral)) * dx(1))
                print(f"J-integral at t={t}: {J_integral}")
                
                # Store results
                fname2.write(str(t*u_r) + "\t")
                fname2.write(str(J_integral) + "\n")

fname2.close()
print('Simulation completed')
