from dolfin import *
from mshr import *
import numpy as np

# Region for J Int
R1_val = 0.25
R2_val = 0.3

# Domain
L = 1
h = 1

# Mesh
domain = Rectangle(Point(-L/2, -h/2), Point(L/2, h/2))
domain.set_subdomain(1, (Circle(Point(0, 0), R2_val) - Circle(Point(0, 0), R1_val) - Rectangle(Point(-R2_val-0.02, -0.01), Point(-R1_val+0.02,0.01))))
domain.set_subdomain(2, Rectangle(Point(-L/2,-0.01), Point(0,0)))
mesh = generate_mesh(domain, 100)
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim(), mesh.domains())
boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)
dx = Measure('dx', domain=mesh, subdomain_data=subdomains)

# Save markers for visualization
File("subdomains.pvd") << subdomains

class TopBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.5)

class BottomBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], -0.5)

def Crack(x):
    return abs(x[1]) < 1e-03 and x[0] <= 0.0

r = Expression(('x[0]','x[1]'),degree=1)


top_boundary = TopBoundary()
bottom_boundary = BottomBoundary()

# Function Spaces
V = FunctionSpace(mesh, 'CG', 1)
W = VectorFunctionSpace(mesh, 'CG', 1)
WW = FunctionSpace(mesh, 'DG', 0)
p, q = TrialFunction(V), TestFunction(V)
FS = FunctionSpace(mesh, 'Lagrange', 1)
P1 = FiniteElement('Lagrange', mesh.ufl_cell(), degree=2)
R1 = FiniteElement('Real', mesh.ufl_cell(), 0)
MFS = FunctionSpace(mesh, MixedElement([(P1*P1),R1,R1]))
f = Function(MFS)
v_f = TestFunction(MFS)
u, c_trans, c_rot = split(f)

# Material parameters
Gc =  1.5
l = 0.015
lmbda = 121.1538e3
mu = 80.7692e3
Temp_expr = Expression("t", t=0.0, degree=1)
# Constitutive functions
def epsilon(u):
    return sym(grad(u)) - as_matrix([[0,0],[0,Temp_expr]])
def sigma(u):
    return 2.0*mu*epsilon(u) + lmbda*tr(epsilon(u))*Identity(len(u))
def psi(u):
    return 0.5 * (lmbda + mu) * (0.5 * (tr(epsilon(u)) + abs(tr(epsilon(u)))))**2 + mu * inner(dev(epsilon(u)), dev(epsilon(u)))		
def H(uold, unew, Hold):
    return conditional(lt(psi(uold), psi(unew)), psi(unew), Hold)
def bulk_energy_density(u):
    return 0.5*inner(sigma(u),epsilon(u))
bc_phi = [DirichletBC(V, Constant(1.0), Crack)]

# Variational form
unew, uold = Function(W), Function(W)
pnew, pold, Hold = Function(V), Function(V), Function(V)

# Strain energy and Phase field energy
E_du = ((1.0 - pold)**2) * bulk_energy_density(u) * dx + (c_trans*u[0])*dx + c_rot*(r[0]*u[1]-r[1]*u[0])*dx
E_phi = (Gc * l * inner(grad(p), grad(q)) + ((Gc / l) + 2.0 * H(uold, unew, Hold)) * inner(p, q) - 2.0 * H(uold, unew, Hold) * q) * dx

Res = derivative(E_du, f, v_f)  # Residual for strain energy


p_phi = LinearVariationalProblem(lhs(E_phi), rhs(E_phi), pnew, bc_phi)


load_top = Constant(0.0)
load_bot = Constant(0.0)
bcbot = DirichletBC(MFS.sub(0).sub(1), load_bot, bottom_boundary)
bctop = DirichletBC(MFS.sub(0).sub(1), load_top, top_boundary)
bc_u = [bcbot, bctop]

# Solver for phase field
solver_phi = LinearVariationalSolver(p_phi)

# Weight function for J Integral Calculation
q_func_space = FunctionSpace(mesh, 'CG', 1)
q_expr = Expression('sqrt(pow(x[0], 2) + pow(x[1], 2)) < r1 ? \
                     1.0 : (sqrt(pow(x[0], 2) + pow(x[1], 2)) > r2 ? 0.0 : \
                     (r2 - sqrt(pow(x[0], 2) + pow(x[1], 2))) / (r2 - r1))', \
                          degree=1, r1=R1_val, r2=R2_val)
q_J_integral = interpolate(q_expr, q_func_space)

# Initialization and simulation loop
t = 0
deltaT = 0.0001  # For load steps
tol = 1e-3        # For staggered scheme
conc_f = File("./ResultsDir/phi.pvd")
conc_u = File("./ResultsDir/u.pvd")
sigmaxx_file = File("./ResultsDir/sxx.pvd")
sigmaxy_file = File("./ResultsDir/sxy.pvd")
sigmayy_file = File("./ResultsDir/syy.pvd")
fname2 = open('JInt.txt', 'w')
fname3 = open('Itercount.txt','w')
while t <= 0.004:
    t += deltaT
    Temp_expr.t = -t
    load_bot.t = -t
    load_top.t = t
    iter, err = 0, 1

    while err > tol:
        iter += 1
        solve(Res == 0, f, bc_u,
              solver_parameters={
                  "newton_solver": {
                      "maximum_iterations": 100,
                      "relative_tolerance": 1e-8,
                      "absolute_tolerance": 1e-10,
                      "linear_solver": "mumps",  # or "lu"
                      "report": True,
                      "error_on_nonconvergence": True,
                      "relaxation_parameter": 0.9
                  }
              })
        u_mixed, _, _ = f.split(deepcopy=True)
        unew.assign(u_mixed)
        solver_phi.solve()
        err_u = errornorm(unew, uold, norm_type='l2')
        err_phi = errornorm(pnew, pold, norm_type='l2')
        err = max(err_u, err_phi)

        uold.assign(unew)
        pold.assign(pnew)
        Hold.assign(project(psi(unew), WW))

        if err < tol:
            print('Iterations:', iter, ', Total time', t)
            sigma_xx = sigma(unew)[0, 0]
            sigma_xx = project(sigma_yy, WW)
            sigma_xx.rename("sigma_xx","sigma xx")
            sigmaxx_file << sigma_xx
            sigma_xy = sigma(unew)[0, 1]
            sigma_xy = project(sigma_xy, WW)
            sigma_xy.rename("sigma_xy","sigma xy")
            sigmaxy_file << sigma_xy
            sigma_yy = sigma(unew)[1, 1]
            sigma_yy = project(sigma_yy, WW)
            sigma_yy.rename("sigma_yy","sigma yy")
            sigmayy_file << sigma_yy
            conc_f << pnew
            conc_u << unew
            fname3.write(str(t)+ "\t")
            fname3.write(str(iter) + "\n")
            J_integral = assemble(dot((sigma(unew)*(grad(unew)) - 0.5*inner(sigma(unew),epsilon(unew))*Identity(2))[:,0] , grad(q_J_integral)) * dx(1))
            print(f"J-integral at t={t}: {J_integral}")
            fname2.write(str(t) + "\t")
            fname2.write(str(J_integral) + "\n")
            
fname3.close()
fname2.close()
print('Simulation completed')

