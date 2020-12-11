from dolfin import (dx, Constant, assemble_system, TestFunction, 
                    as_backend_type, assemble, PETScOptions, Function, plot, File, 
                    PETScMatrix, PETScVector, MPI)
from matplotlib.pyplot import subplots
import petsc4py
from slepc4py import SLEPc
from petsc4py import PETSc
import numpy as np
import dolfin
import ufl

class EigenSolver(object):
    def __init__(self, 
                 a_k, 
                 u,
                 a_m=None,                                   # optional, for eigpb of the type (K-lambda M)x=0
                 bcs=None,
                 restricted_dofs_is=None,
                 slepc_options={'eps_max_it':100},
                 option_prefix=None,
                 comm=MPI.comm_world, 
                 slepc_eigensolver = None
                ):
        self.comm = comm
        self.slepc_options = slepc_options
        if option_prefix:
            self.E.setOptionsPrefix(option_prefix)
        self.V = u.function_space()
        self.index_set_not_bc = None
        if type(bcs) == list:
            self.bcs = bcs
        elif type(bcs) == dolfin.fem.dirichletbc.DirichletBC:
            self.bcs = [bcs]
        else:
            self.bcs = None

        if type(a_k) == ufl.form.Form:
            # a form to be assembled
            self.K = as_backend_type(assemble(a_k)).mat()
        elif type(a_k) == petsc4py.PETSc.Mat:
            # an assembled petsc matrix
            self.K = a_k

        if a_m is not None and type(a_m) == ufl.form.Form:
            self.M = as_backend_type(assemble(a_m)).mat()
        elif a_m is not None and type(a_m) == petsc4py.PETSc.Mat:
            self.M = a_m

        # if bcs extract reduced matrices on dofs with no bcs
        if self.bcs:
            self.index_set_not_bc = self.get_interior_index_set(self.bcs, self.V)
        elif restricted_dofs_is:
            self.index_set_not_bc = restricted_dofs_is

        if self.index_set_not_bc is not None:
            try:
                self.K = self.K.createSubMatrix(self.index_set_not_bc, self.index_set_not_bc)
                if a_m is not None:
                   self.M = self.M.createSubMatrix(self.index_set_not_bc, self.index_set_not_bc)
            except:
                self.K = self.K.getSubMatrix(self.index_set_not_bc, self.index_set_not_bc)
                if a_m is not None:
                    self.M = self.M.getSubMatrix(self.index_set_not_bc, self.index_set_not_bc)
            self.projector = petsc4py.PETSc.Scatter()
            self.projector.create(
                vec_from=self.K.createVecRight(),
                is_from=None,
                vec_to=u.vector().vec(),
                is_to=self.index_set_not_bc
                )

        # set up the eigensolver
        if slepc_eigensolver:
            self.E = slepc_eigensolver
        else:
            self.E = self.eigensolver_setup()

        if a_m is not None:
            print(a_m)
            self.E.setOperators(self.K,self.M)       
        else:
            self.E.setOperators(self.K)
    
    def get_interior_index_set(self, boundary_conditions, function_space):
        """Returns the index set with free dofs"""
        # Find dofs affected by the boundary conditions
        bc_dofs = []    
        for bc in boundary_conditions:
            bc_dofs.extend(bc.get_boundary_values().keys()) 
        ownership_range = function_space.dofmap().ownership_range()
        interior_dofs = [x for x in range(ownership_range[0], ownership_range[1]) if x not in bc_dofs]    
        # Create petsc4py.PETSc.IS object with interior degrees of freedom
        index_set = petsc4py.PETSc.IS()
        index_set.createGeneral(interior_dofs)  
        return index_set

    def eigensolver_setup(self):
        E = SLEPc.EPS()
        E.create()        
        E.setType(SLEPc.EPS.Type.KRYLOVSCHUR)
        #E.setType(SLEPc.EPS.Type.LANCZOS)
        E.setProblemType(SLEPc.EPS.ProblemType.GHIEP)
        #E.setWhichEigenpairs(E.Which.LARGEST_MAGNITUDE)
        E.setWhichEigenpairs(E.Which.TARGET_MAGNITUDE)
        E.setTarget(0.) 
        st = E.getST()
        st.setType('sinvert')
        st.setShift(1.e-3)
        return E
        
    
    def solve(self, n_eig, verbose=True):
        E = self.E
        E.setDimensions(n_eig)
        self.set_options(self.slepc_options,verbose=verbose)
        E.setFromOptions()
        E.solve()
        # print info
        its = E.getIterationNumber()
        eps_type = E.getType()
        self.nev, ncv, mpd = E.getDimensions()
        tol, maxit = E.getTolerances()
        self.nconv = E.getConverged()
        if verbose == True:
            print("Solution method: {:s}, stopping condition: tol={:.4g}, maxit={:d}".format(eps_type,tol, maxit))
            print("Number of converged/requested eigenvalues with {:d} iterations  : {:d}/{:d}".format(its,self.nconv,self.nev))
        return self.nconv, its
        
    def set_options(self,slepc_options,verbose=True):
        if verbose == True:
            print("---- setting additional slepc options -----")
        for (opt, value) in slepc_options.items():
            if verbose == True:
                print("    ",opt,":",value)
            PETScOptions.set(opt,value) 
        self.E.setFromOptions()
        
    def get_eigenpair(self,i):
        u_r = Function(self.V)
        u_im = Function(self.V)
        v_r, v_i = self.K.createVecs()
        eig = self.E.getEigenpair(i, v_r, v_i)
        err = self.E.computeError(i)
        if self.index_set_not_bc:
            self.projector.scatter(vec_from=v_r, vec_to=u_r.vector().vec())
            self.projector.scatter(vec_from=v_i, vec_to=u_im.vector().vec())
            u_r.vector().vec().ghostUpdate()
            u_im.vector().vec().ghostUpdate()
        return eig, u_r, u_im, err
        
    def get_eigenvalues(self,n):
        eigenvalues = [] 
        for i in range(n):
            eig, u_r, u_im, err = self.get_eigenpair(i)
            eigenvalues.append([eig.real, err])
        return np.array(eigenvalues)
    
    def get_eigenpairs(self,n):
        eigenvalues = [] 
        eigenvectors_real = [] 
        eigenvectors_im = [] 
        for i in range(n):
            eig, u_r, u_im, err = self.get_eigenpair(i)
            eigenvalues.append(eig)
            eigenvectors_real.append(u_r)
            eigenvectors_im.append(u_im)
        return np.array(eigenvalues), [eigenvectors_real, eigenvectors_im]
    
    def save_eigenvectors(self,n,file_name="output/modes.pvd",save_imaginary=False):
        eigenvalues = [] 
        eigenvectors = [] 
        file = File(self.comm,file_name)
        for i in range(n):
            eig, u_r, u_im, err = self.get_eigenpair(i)
            u_r.rename("mode real","mode real")
            file.write(u_r,i)
            if save_imaginary:
                u_im.rename("mode imaginary","mode imaginary")
                file.write(u_im,i)
        print('saved eigenvectors in ', file_name)
        return file_name
    
    def plot_eigenpair(self,i):
        eig, u_r, u_im, err = self.get_eigenpair(i)
        p = plot(u_r,title="{:d} -- {:2.2f}".format(i,eig))
        return p 
    

if __name__ == "__main__":

    import dolfin
    import ufl
    import numpy as np
    import matplotlib.pyplot as plt
    dolfin.parameters["use_petsc_signal_handler"] = True
    dolfin.parameters["form_compiler"]["cpp_optimize"] = True
    dolfin.parameters["form_compiler"]["representation"] = "uflacs"

    n = 100
    mesh = dolfin.UnitSquareMesh(n,n)
    V = dolfin.FunctionSpace(mesh,'CG',1)
    u = dolfin.Function(V)
    ut = dolfin.TestFunction(V)
    v = dolfin.TrialFunction(V)
    dx = dolfin.Measure("dx",domain=mesh)
    bc = dolfin.DirichletBC(V,0,"on_boundary") 
    a_k = ufl.dot(ufl.grad(ut),ufl.grad(v))*dx
    a_m = ut*v*dx
    eig_solver = EigenSolver(a_k, u, a_m, [bc])
    #eig_solver = EigenSolver(a_k, u, bcs=[bc])
    plt.savefig("operators.png")
    ncv, it = eig_solver.solve(10)
    eigs = eig_solver.get_eigenvalues(ncv)
    plt.figure()
    plt.plot(eigs,'o')
    plt.title("Eigenvalues")
    plt.savefig("eigenvalues.png")

    eig_solver.save_eigenvectors(ncv)
    for i in range(ncv): 
        plt.figure()
        eig_solver.plot_eigenpair(i)
        plt.savefig("output/mode-{:d}.png".format(i))


