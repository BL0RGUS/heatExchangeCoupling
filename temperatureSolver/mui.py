import mui4py
from mpi4py import MPI
import numpy as np


class MUI:

    def __init__(self, nodes, dt):
        self.nodes = nodes
        self.dt = dt

        # Interface setup
        mui4py.mpi_split_by_app()

        self.MPI_COMM_WORLD = MPI.COMM_WORLD

        dims = 1
        config = mui4py.Config(dims, mui4py.FLOAT64)
        
        self.solverNum = self.MPI_COMM_WORLD.Get_rank()
        self.numSolvers = self.MPI_COMM_WORLD.Get_size()


        if self.numSolvers == 1:
            return None

        if self.solverNum == 0:
            iface =  ["ifs1"]
        elif self.solverNum == self.numSolvers-1:
            iface =  ["ifs" + str(self.solverNum)]
        else:
            iface = ["ifs"+str(self.solverNum), "ifs"+str(self.solverNum+1)]

        domain = "Solver" + str(self.solverNum)
        unifaces = mui4py.create_unifaces(domain, iface, config)

        if self.solverNum != 0:
            self.leftUniface = unifaces["ifs" + str(self.solverNum)]
            self.leftUniface.set_data_types({"temp": mui4py.FLOAT64,
                                     "alpha": mui4py.FLOAT64})

        if self.solverNum != self.numSolvers-1:
            self.rightUniface = unifaces["ifs" + str(self.solverNum+1)]

            self.rightUniface.set_data_types({"temp": mui4py.FLOAT64,
                                     "alpha": mui4py.FLOAT64,
                                     "nodes":mui4py.INT})

        # Use MPI allreduce to find the solver with the smallest dt and the most nodes
        self.minDt, self.maxNodes = self.findSuperlativeParameters()

        #TODO: change spatial sampler so we can use different mesh resolution between solvers
        self.s_sampler = mui4py.SamplerPseudoNearestNeighbor(0.5)
        self.t_sampler = mui4py.TemporalSamplerExact()

        ## Point array for fetching using fetch_many and push_many
        # needs to be normalised to the maximum number nodes for any solver
        i = 0
        c = 0
        step = self.maxNodes/self.nodes
        self.points = np.zeros((self.nodes, 1))
        while c < self.nodes:
            self.points[c] = [i]
            i += step
            c += 1

    def findSuperlativeParameters(self):
        minDt = self.MPI_COMM_WORLD.allreduce(self.dt, op=MPI.MIN)
        maxNodes = self.MPI_COMM_WORLD.allreduce(self.nodes, op=MPI.MAX)
        
        return (minDt, maxNodes)

    def getAlphas(self, alpha):
        rightAlpha = None
        leftAlpha = None
        if self.solverNum != 0:
            self.leftUniface.push("alpha", [self.solverNum], alpha)
            self.leftUniface.commit( 0 )
            leftAlpha = self.leftUniface.fetch("alpha", [self.solverNum-1], 0, mui4py.SamplerExact(), mui4py.TemporalSamplerExact())
        else:
            leftAlpha = alpha
        if self.solverNum != self.numSolvers-1:
            self.rightUniface.push("alpha", [self.solverNum], alpha)
            self.rightUniface.commit( 0 )

            rightAlpha = self.rightUniface.fetch("alpha", [self.solverNum+1], 0, mui4py.SamplerExact(), mui4py.TemporalSamplerExact())
        else:
            rightAlpha = alpha
            
        return (leftAlpha, rightAlpha)

    
    def pushRight(self, vals, time, data="temp"):
        self.rightUniface.push_many(data, self.points, vals)
        self.rightUniface.commit( time )

    def pushLeft(self, vals, time, data="temp"):
        self.leftUniface.push_many(data, self.points, vals)
        self.leftUniface.commit( time )
            

    def fetchRightPrev(self, time, data="temp"):
        vals = self.leftUniface.fetch_many(data, self.points, time, 
                                       self.s_sampler, self.t_sampler)
        
        # forget to save memory
        self.leftUniface.forget( time-self.dt)
        #print( vals )
        return vals
    
    def fetchLeftNext(self, time, data="temp"):
        vals = self.rightUniface.fetch_many(data, self.points, time, 
                                       self.s_sampler, self.t_sampler)
        
        # forget to save memory
        self.rightUniface.forget( time-self.dt)
        #print( vals )
        return vals
